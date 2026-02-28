package iso;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.NavigableMap;
import java.util.TreeMap;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.ejml.data.DMatrixRMaj;
import org.jgrapht.Graph;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.alg.connectivity.ConnectivityInspector;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.StandardOpenOption;
import java.util.stream.Stream;
import java.util.concurrent.atomic.AtomicLong;

public class Main
{
    public static int COEFFICIENTS = 6;
    public static int N_THREADS = 4;
    public static int ORIG_GADGET_SIZE = 3;
    public static int ORIG_SUBD = 1;
    public static int MAX_ORIG_NODES = 500;
    public static int MIN_VERTICES_AFTER_SUBDIVISION = 50;
    public static int MAX_VERTICES_AFTER_SUBDIVISION = 1500;
    public static long SEED = 12345L;
    public static boolean USE_TRIPLETS_ON_PAIR_FAILURE = true;
    public static boolean SKIP_TRIPLETS_ON_FIRST_ROUND = true;
    public static int CURRENT_ROUND = 0;

    public static boolean DEBUG = true;
    public static int TOP_OUTCOMES_TO_SHOW = 1;
    /** Generic progress log time gate (ms). */
    public static long PROGRESS_LOG_EVERY_MS = 2000; // adjust
    public static AtomicLong LAST_PROGRESS_LOG_MS = new AtomicLong(0);
    /** Global probe progress log time gate (ms). */
    public static long PROBE_PROGRESS_LOG_EVERY_MS = 2000; // adjust
    public static AtomicLong LAST_PROBE_PROGRESS_LOG_MS = new AtomicLong(0);

    public static String GRAPH_PATH = "graph.txt";
    public static String GRAPH_FOLDER = "/graph_folder";
    public static String OUTPUT_LOG = "/graph_folder/OUTPUT_LOGS_SEED_12345.csv";
    public static long MAX_INPUT_BYTES = 100L * 1024L * 1024L;
    public static int MAX_ORIG_EDGES = 10_000;

    public static double EIGEN_DEFAULT_SCALE_OV = 1e6; //1e10
    public static double DS_EPS_OV = 1e-6; //1e-10;
    public static double CURV_COEFF_DEFAULT_SCALE_OV = 1e6;  //1e10
    
    public enum TwinKind { NONE, TRUE_TWINS, FALSE_TWINS }

    public static void main(String[] args)
    {
    	applyCliOverrides(args);
        Path root = Path.of(GRAPH_FOLDER);
        Path baseLog = Path.of(OUTPUT_LOG);
        Path successLog = deriveLogPath(baseLog, "SUCCESS");
        Path failLog = deriveLogPath(baseLog, "FAIL");
        Path exceptionLog = deriveLogPath(baseLog, "EXCEPTION");

        try {
            Path parent = baseLog.getParent();
            if (parent != null) Files.createDirectories(parent);
            ensureCsvHeader(successLog);
            ensureCsvHeader(failLog);
            ensureCsvHeader(exceptionLog);
        } catch (IOException e) {
            e.printStackTrace();
            return;
        }

        final java.util.Set<String> processed = new java.util.HashSet<>();
        try {
            processed.addAll(readFirstColumnCsv(successLog));
            processed.addAll(readFirstColumnCsv(failLog));
            processed.addAll(readFirstColumnCsv(exceptionLog));
            // remove header token if present
            processed.remove("graph_path");
        } catch (IOException e) {
            e.printStackTrace();
            return;
        }

        try (Stream<Path> paths = Files.walk(root)) {
            paths
                .filter(Files::isRegularFile)
                .filter(p -> !p.equals(successLog) && !p.equals(failLog) && !p.equals(exceptionLog))
                .forEach(p -> {
                    String graphPath = p.toString();

                    if (processed.contains(graphPath)) {
                        if (DEBUG) System.out.println("Skipping already-processed: " + graphPath);
                        return;
                    }

                    try {
                        boolean solved = solve(graphPath);
                        if (solved) {
                            System.out.println(graphPath + " SOLVED SUCCESSFULLY.");
                            appendCsvLineQuiet(successLog, graphPath, "SUCCESS", "");
                        } else {
                            System.out.println(graphPath + " FAILED.");
                            appendCsvLineQuiet(failLog, graphPath, "FAIL", "returned false");
                        }
                        processed.add(graphPath); // so duplicates in folder walk are skipped too
                    } catch (Exception ex) {
                        ex.printStackTrace();
                        appendCsvLineQuiet(exceptionLog, graphPath, "EXCEPTION", safeMessage(ex));
                        processed.add(graphPath);
                    }
                });
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static boolean solve(String graphPath) throws Exception
    {
    	Constants.EIGEN_DEFAULT_SCALE = EIGEN_DEFAULT_SCALE_OV;
    	Constants.DS_EPS = DS_EPS_OV;
    	Constants.CURV_COEFF_DEFAULT_SCALE = CURV_COEFF_DEFAULT_SCALE_OV;  //1e10
        
        Path graphFile = Path.of(graphPath);

        if (isTooHugeToLoad(graphFile)) {
            throw new IllegalStateException("Too huge to load (size/header limits)");
        }

        Graph<Integer, DefaultEdge> g1raw = GraphUtils.readGraphAuto(graphFile);
        final int nOriginal = g1raw.vertexSet().size();

        if (nOriginal > MAX_ORIG_NODES) {
            throw new IllegalStateException(
                "Original graph too large. " +
                "Limit=" + MAX_ORIG_NODES +
                "  nOriginal=" + nOriginal +
                "  graphPath=" + graphPath
            );
        }

        if (!new ConnectivityInspector<>(g1raw).isConnected()) {
            throw new IllegalArgumentException("Input graph is not connected." + graphPath);
        }

        int[] perm = GraphUtils.randomPermutation(nOriginal, SEED);
        Graph<Integer, DefaultEdge> g2raw = GraphUtils.permute(g1raw, perm);

        // Separate working copies for augmentation:
        Graph<Integer, DefaultEdge> g1base = GraphUtils.readGraphAuto(graphFile);
        Graph<Integer, DefaultEdge> g2base = GraphUtils.permute(g1base, perm);

        System.out.println("Graph: " + graphPath);
        System.out.println("Original nodes: " + nOriginal);

        int subd = ORIG_SUBD;
        Graph<Integer, DefaultEdge> g1 = g1base;
        Graph<Integer, DefaultEdge> g2 = g2base;

        while (true) {
            Graph<Integer, DefaultEdge> g1try = GraphUtils.subdivideEveryEdge(g1base, subd);
            Graph<Integer, DefaultEdge> g2try = GraphUtils.subdivideEveryEdge(g2base, subd);

            int g1SizeTry = g1try.vertexSet().size();
            int g2SizeTry = g2try.vertexSet().size();

            if (g1SizeTry > MAX_VERTICES_AFTER_SUBDIVISION || g2SizeTry > MAX_VERTICES_AFTER_SUBDIVISION) {
                throw new IllegalStateException(
                        "Graph too large after subdivision. " +
                        "Limit=" + MAX_VERTICES_AFTER_SUBDIVISION +
                        "  subd=" + subd +
                        "  g1Vertices=" + g1SizeTry +
                        "  g2Vertices=" + g2SizeTry
                );
            }

            boolean bigEnough = (g1SizeTry >= MIN_VERTICES_AFTER_SUBDIVISION && g2SizeTry >= MIN_VERTICES_AFTER_SUBDIVISION);

            // Only bother testing Ds once we’re big enough (optional but usually faster).
            boolean dsOk = false;
            if (bigEnough) {
                boolean ds1ok = spectralDimensionSucceeds(g1try);
                boolean ds2ok = spectralDimensionSucceeds(g2try);
                dsOk = ds1ok && ds2ok;

                if (DEBUG) {
                    System.out.println("Subdivision trial subd=" + subd
                            + "  sizes=(" + g1SizeTry + "," + g2SizeTry + ")"
                            + "  DsOk=(" + ds1ok + "," + ds2ok + ")");
                }
            } else {
                if (DEBUG) {
                    System.out.println("Subdivision trial subd=" + subd
                            + "  sizes=(" + g1SizeTry + "," + g2SizeTry + ")"
                            + "  below MIN=" + MIN_VERTICES_AFTER_SUBDIVISION);
                }
            }

            if (bigEnough && dsOk) {
                g1 = g1try;
                g2 = g2try;
                System.out.println("Subdivision chosen: subd=" + subd + " -> vertices=" + g1SizeTry + " edges= " + g1.edgeSet().size()
                        + " (spectral dimension succeeded)");
                break;
            }

            // Not acceptable yet -> subdivide more
            subd++;
        }

        int g1Size = g1.vertexSet().size();
        int g2Size = g2.vertexSet().size();

        if (g1Size > MAX_VERTICES_AFTER_SUBDIVISION || g2Size > MAX_VERTICES_AFTER_SUBDIVISION) {
            throw new IllegalStateException(
                "Graph too large after subdivision. " +
                "Limit=" + MAX_VERTICES_AFTER_SUBDIVISION +
                "  g1Vertices=" + g1Size +
                "  g2Vertices=" + g2Size
            );
        }

        // Start gadget size at configured ORIG_GADGET_SIZE and increment after each *matched pair* attachment.
        int gadgetSize = ORIG_GADGET_SIZE;

        // --- helper: build reverse mapping G2->G1 from groupIndexG1ToG2 ---
        java.util.function.Function<TwoGraphPartitionResult, int[]> buildReverseMapG2toG1 = (refined) -> {
            int[] mapG2toG1 = new int[refined.groupsG2.size()];
            Arrays.fill(mapG2toG1, -1);
            for (int i = 0; i < refined.groupIndexG1ToG2.length; i++) {
                int j = refined.groupIndexG1ToG2[i];
                if (j >= 0) mapG2toG1[j] = i;
            }
            return mapG2toG1;
        };

        // --- helper: collect all mapped singleton pairs (u in G1, v in G2), dedup + deterministic sort ---
        java.util.function.Function<TwoGraphPartitionResult, List<int[]>> collectMappedSingletonPairs = (refined) -> {
            int[] mapG2toG1 = buildReverseMapG2toG1.apply(refined);
            List<int[]> pairs = new ArrayList<>();

            // From G1 -> G2
            for (int i = 0; i < refined.groupsG1.size(); i++) {
                if (refined.groupsG1.get(i).size() != 1) continue;
                int j = refined.groupIndexG1ToG2[i];
                if (j < 0) continue;
                if (refined.groupsG2.get(j).size() != 1) continue;

                int u = refined.groupsG1.get(i).get(0);
                int v = refined.groupsG2.get(j).get(0);
                pairs.add(new int[]{u, v});
            }

            // From G2 -> G1 (covers the other direction too)
            for (int j = 0; j < refined.groupsG2.size(); j++) {
                if (refined.groupsG2.get(j).size() != 1) continue;
                int i = mapG2toG1[j];
                if (i < 0) continue;
                if (refined.groupsG1.get(i).size() != 1) continue;

                int u = refined.groupsG1.get(i).get(0);
                int v = refined.groupsG2.get(j).get(0);
                pairs.add(new int[]{u, v});
            }

            // Deterministic sort + dedup
            pairs.sort((a, b) -> {
                int c = Integer.compare(a[0], b[0]);
                if (c != 0) return c;
                return Integer.compare(a[1], b[1]);
            });

            List<int[]> dedup = new ArrayList<>();
            int prevU = Integer.MIN_VALUE, prevV = Integer.MIN_VALUE;
            for (int[] p : pairs) {
                if (p[0] != prevU || p[1] != prevV) dedup.add(p);
                prevU = p[0]; prevV = p[1];
            }
            return dedup;
        };

        // --- helper: pick smallest *mapped* refined group pair (i in G1, j in G2) deterministically ---
        java.util.function.Function<TwoGraphPartitionResult, int[]> pickSmallestMappedGroupPair = (refined) -> {
            int bestI = -1;
            int bestJ = -1;
            int bestScore = Integer.MAX_VALUE;

            for (int i = 0; i < refined.groupsG1.size(); i++) {
                int j = refined.groupIndexG1ToG2[i];
                if (j < 0) continue;

                int sz1 = refined.groupsG1.get(i).size();
                int sz2 = refined.groupsG2.get(j).size();
                int score = Math.min(sz1, sz2);

                if (score < bestScore) {
                    bestScore = score;
                    bestI = i;
                    bestJ = j;
                } else if (score == bestScore && bestI >= 0) {
                    // Tie-breaker: ProbeProfile ordering on the G1 side (deterministic)
                    ProbeProfile pA = refined.groupProfilesG1.get(i);
                    ProbeProfile pB = refined.groupProfilesG1.get(bestI);
                    if (pA.compareTo(pB) < 0) {
                        bestI = i;
                        bestJ = j;
                    }
                }
            }

            if (bestI < 0) {
                throw new IllegalStateException("No mapped refined groups found across graphs (all G1 groups unmapped).");
            }
            return new int[]{bestI, bestJ};
        };

        int round = 0;
        while (true)
        {
            round++;
            CURRENT_ROUND = round;
            System.out.println("========== ROUND " + round + " ==========");

            // --- Compute curvature coefficients for both graphs ---
            CurvatureUtils.Result c1;
            CurvatureUtils.Result c2;

            {
                int n1 = g1.vertexSet().size();
                DMatrixRMaj V1 = new DMatrixRMaj(n1, n1);
                double[] evals1 = EigenUtils.computeSpectrumFromGraph(g1, EigenUtils.GraphMatrixType.LAPLACIAN, V1);
                SpectralUtils.DsResult ds1 = SpectralUtils.estimateDs(evals1, Constants.SPECTRAL_EXPONENT, Constants.DIFFUSION_CONSTANT);
                c1 = CurvatureUtils.compute(evals1, V1, ds1.ds, COEFFICIENTS);
            }

            {
                int n2 = g2.vertexSet().size();
                DMatrixRMaj V2 = new DMatrixRMaj(n2, n2);
                double[] evals2 = EigenUtils.computeSpectrumFromGraph(g2, EigenUtils.GraphMatrixType.LAPLACIAN, V2);
                SpectralUtils.DsResult ds2 = SpectralUtils.estimateDs(evals2, Constants.SPECTRAL_EXPONENT, Constants.DIFFUSION_CONSTANT);
                c2 = CurvatureUtils.compute(evals2, V2, ds2.ds, COEFFICIENTS);
            }

            boolean stopBySingletonOrTwins =
                    allClassesSingletonOrTwinsAndMatchAcrossGraphs(
                            g1, c1.qCoeffs,
                            g2, c2.qCoeffs,
                            g1raw, g2raw,
                            nOriginal,
                            Constants.SIGNATURE_BFS_LAYERS
                    );

            System.out.println("Stop condition (singletons or twin-classes)? " + stopBySingletonOrTwins);

            if (stopBySingletonOrTwins) {
                System.out.println("DONE: all original classes are singletons or twin-classes (matching across graphs).");

                int[] map = buildOriginalMappingByTypeAllowingTwins(
                        g1, c1.qCoeffs,
                        g2, c2.qCoeffs,
                        g1raw, g2raw,
                        nOriginal,
                        Constants.SIGNATURE_BFS_LAYERS
                );

                boolean ok = isIsomorphismOnRawGraphs(g1raw, g2raw, map, nOriginal);
                System.out.println("MAPPING-ISOMORPHISM check on RAW original graphs: " + ok);
                return ok;
            }

            // Pick the smallest common original-only equivalence class (size >= 2) across the two graphs.
            List<List<Integer>> smallest =
                    smallestCommonOriginalEquivalenceClassVertsAtLeast2(
                            g1, c1.qCoeffs,
                            g2, c2.qCoeffs,
                            nOriginal,
                            Constants.SIGNATURE_BFS_LAYERS
                    );

            if (smallest == null)
                throw new IllegalStateException("No matching non-singleton class found across graphs (>=2).");

            List<Integer> candidates1 = smallest.get(0);
            List<Integer> candidates2 = smallest.get(1);

            System.out.println("Smallest matching class candidates:");
            System.out.println("  G1: " + candidates1);
            System.out.println("  G2: " + candidates2);

            // Marker paths (must differ)
            int pathLenU = gadgetSize + 1;
            int pathLenV = gadgetSize + 2;

            TwoGraphPartitionResult refined =
                    partitionVerticesBySharedCliqueProbingTwoGraphs(
                            g1, candidates1, nOriginal,
                            g2, candidates2, nOriginal,
                            /*K*/ gadgetSize,
                            pathLenU, pathLenV,
                            DEBUG,
                            TOP_OUTCOMES_TO_SHOW
                    );

            boolean refinedHasAnySingleton1 = refined.groupsG1.stream().anyMatch(gr -> gr.size() == 1);
            boolean refinedHasAnySingleton2 = refined.groupsG2.stream().anyMatch(gr -> gr.size() == 1);

            boolean refinedNoSingletons1 = refined.groupsG1.stream().noneMatch(gr -> gr.size() == 1);
            boolean refinedNoSingletons2 = refined.groupsG2.stream().noneMatch(gr -> gr.size() == 1);

            // Case A: NO singletons on both sides => pick from the SAME mapped refined group pair.
            if (refinedNoSingletons1 && refinedNoSingletons2)
            {
                int[] ij = pickSmallestMappedGroupPair.apply(refined);
                int i = ij[0];
                int j = ij[1];

                List<Integer> g1Group = refined.groupsG1.get(i);
                List<Integer> g2Group = refined.groupsG2.get(j);

                int pick1 = g1Group.get(0);
                int pick2 = g2Group.get(0);

                System.out.println("Refinement returned NO singletons on both sides.");
                System.out.println("Picked vertices from the SAME mapped refined group:");
                System.out.println("  G1 groupIndex=" + i + " pick=" + pick1 + " from group=" + g1Group);
                System.out.println("  G2 groupIndex=" + j + " pick=" + pick2 + " from group=" + g2Group);

                g1 = GraphUtils.attachCliqueAsSuperClique(g1, pick1, gadgetSize);
                g2 = GraphUtils.attachCliqueAsSuperClique(g2, pick2, gadgetSize);
                if (DEBUG) {
                    System.out.println("Attached superclique(size=" + gadgetSize + ") to mapped picks: G1 " + pick1 + " <-> G2 " + pick2);
                }

                gadgetSize++;
                System.out.println();
                continue;
            }

            // Case B: There exist singleton(s) somewhere => attach ONLY mapped singleton pairs (same refined group across graphs).
            if (refinedHasAnySingleton1 || refinedHasAnySingleton2)
            {
                List<int[]> singletonPairs = collectMappedSingletonPairs.apply(refined);

                if (singletonPairs.isEmpty()) {
                    // This means you have singletons but none that correspond across graphs (mapping mismatch).
                    throw new IllegalStateException(
                            "Singleton(s) exist after refinement, but NO mapped singleton pairs were found across graphs."
                    );
                }

                System.out.println("Refinement found singleton(s). Attaching cliques to ALL mapped singleton pairs found...");

                int attachedPairs = 0;
                for (int[] p : singletonPairs) {
                    int u = p[0];
                    int v = p[1];

                    g1 = GraphUtils.attachCliqueAsSuperClique(g1, u, gadgetSize);
                    g2 = GraphUtils.attachCliqueAsSuperClique(g2, v, gadgetSize);

                    if (DEBUG) {
                        System.out.println("  Attached superclique(size=" + gadgetSize + ") to mapped singleton pair: G1 " + u + " <-> G2 " + v);
                    }

                    gadgetSize++;
                    attachedPairs++;
                }

                System.out.println("Attached " + attachedPairs + " mapped singleton pair clique(s).");
                System.out.println();
                continue;
            }

            // Defensive guard.
            throw new IllegalStateException("Unexpected refinement state: neither no-singletons nor has-singletons.");
        }
    }

    /** Holds the original-only signature partition and its derived structures. */
    private static final class PartitionInfo {
        final NavigableMap<SignatureUtils.VertexSignature, List<Integer>> sigToVertsOriginalOnly;
        final SignatureUtils.VertexSignature[] originalVertexTypesSortedMultiset; // strengthened (multiset of types)
        final int[] originalClassSizesSortedMultiset;                             // strengthened (multiset of class sizes)
        final long[] quantizedSpectrumSorted;                                     // FIX: long[] (probed graph spectral summary)

        PartitionInfo(
                NavigableMap<SignatureUtils.VertexSignature, List<Integer>> sigToVertsOriginalOnly,
                SignatureUtils.VertexSignature[] originalVertexTypesSortedMultiset,
                int[] originalClassSizesSortedMultiset,
                long[] quantizedSpectrumSorted
        ) {
            this.sigToVertsOriginalOnly = sigToVertsOriginalOnly;
            this.originalVertexTypesSortedMultiset = originalVertexTypesSortedMultiset;
            this.originalClassSizesSortedMultiset = originalClassSizesSortedMultiset;
            this.quantizedSpectrumSorted = quantizedSpectrumSorted;
        }
    }

    private static PartitionInfo computePartitionInfoOriginalOnly(Graph<Integer, DefaultEdge> g, int nOriginal)
    {
        int n = g.vertexSet().size();

        DMatrixRMaj V = new DMatrixRMaj(n, n);
        double[] evals = EigenUtils.computeSpectrumFromGraph(g, EigenUtils.GraphMatrixType.LAPLACIAN, V);

        SpectralUtils.DsResult ds = SpectralUtils.estimateDs(
                evals, Constants.SPECTRAL_EXPONENT, Constants.DIFFUSION_CONSTANT
        );

        CurvatureUtils.Result curv = CurvatureUtils.compute(evals, V, ds.ds, COEFFICIENTS);

        NavigableMap<SignatureUtils.VertexSignature, List<Integer>> sigToVerts =
                SignatureUtils.signatureToVerticesOriginalOnly(
                        g, curv.qCoeffs, nOriginal, Constants.SIGNATURE_BFS_LAYERS
                );

        // Multiset of original vertex TYPES (sorted)
        List<SignatureUtils.VertexSignature> typeList = new ArrayList<>(nOriginal);
        // Multiset of original class sizes (sorted)
        List<Integer> sizeList = new ArrayList<>(sigToVerts.size());

        for (var e : sigToVerts.entrySet()) {
            SignatureUtils.VertexSignature type = e.getKey();
            List<Integer> verts = e.getValue();
            sizeList.add(verts.size());
            for (int i = 0; i < verts.size(); i++) typeList.add(type);
        }

        SignatureUtils.VertexSignature[] typesSorted = typeList.toArray(new SignatureUtils.VertexSignature[0]);
        Arrays.sort(typesSorted, SignatureUtils.VertexSignature::compareTo);

        int[] classSizesSorted = new int[sizeList.size()];
        for (int i = 0; i < sizeList.size(); i++) classSizesSorted[i] = sizeList.get(i);
        Arrays.sort(classSizesSorted);

        // Quantized spectrum summary (sorted)  // FIX: long[]
        long[] qspec = EigenUtils.quantize(evals);
        Arrays.sort(qspec);

        return new PartitionInfo(sigToVerts, typesSorted, classSizesSorted, qspec);
    }

    /** Find the signature TYPE that contains original vertex v, plus the class size. */
    private static SigAndSize findSignatureAndSizeOfVertex(
            NavigableMap<SignatureUtils.VertexSignature, List<Integer>> sigToVerts,
            int v
    ) {
        for (var e : sigToVerts.entrySet()) {
            if (e.getValue().contains(v)) {
                return new SigAndSize(e.getKey(), e.getValue().size());
            }
        }
        throw new IllegalStateException("Vertex " + v + " not found in original-only signature partition");
    }

    private static final class SigAndSize {
        final SignatureUtils.VertexSignature sig;
        final int size;
        SigAndSize(SignatureUtils.VertexSignature sig, int size) {
            this.sig = sig;
            this.size = size;
        }
    }

    /** Pair gadget. */
    private static Graph<Integer, DefaultEdge> buildStrongProbeGraphPair(
            Graph<Integer, DefaultEdge> base,
            int u,
            int v,
            int K,
            int uPrivateClique,
            int vPrivateClique,
            int uPathLen,
            int vPathLen
    ) {
        Graph<Integer, DefaultEdge> h = base;

        h = GraphUtils.attachSharedClique(h, Arrays.asList(u, v), K);
        h = GraphUtils.attachCliqueAsSuperClique(h, u, uPrivateClique);
        h = GraphUtils.attachCliqueAsSuperClique(h, v, vPrivateClique);

        if (uPathLen > 0) h = GraphUtils.attachPath(h, u, uPathLen);
        if (vPathLen > 0) h = GraphUtils.attachPath(h, v, vPathLen);

        return h;
    }

    /** Triplet gadget. */
    private static Graph<Integer, DefaultEdge> buildStrongProbeGraphTriplet(
            Graph<Integer, DefaultEdge> base,
            int u,
            int v,
            int w,
            int K,
            int uPrivateClique,
            int vPrivateClique,
            int wPrivateClique,
            int uPathLen,
            int vPathLen,
            int wPathLen
    ) {
        Graph<Integer, DefaultEdge> h = base;

        h = GraphUtils.attachSharedClique(h, Arrays.asList(u, v, w), K);
        h = GraphUtils.attachCliqueAsSuperClique(h, u, uPrivateClique);
        h = GraphUtils.attachCliqueAsSuperClique(h, v, vPrivateClique);
        h = GraphUtils.attachCliqueAsSuperClique(h, w, wPrivateClique);

        if (uPathLen > 0) h = GraphUtils.attachPath(h, u, uPathLen);
        if (vPathLen > 0) h = GraphUtils.attachPath(h, v, vPathLen);
        if (wPathLen > 0) h = GraphUtils.attachPath(h, w, wPathLen);

        return h;
    }

    private static final class EndDesc implements Comparable<EndDesc> {
        final SignatureUtils.VertexSignature type;
        final int size;
        EndDesc(SignatureUtils.VertexSignature t, int s) { this.type = t; this.size = s; }

        @Override public int compareTo(EndDesc o) {
            int c = this.type.compareTo(o.type);
            if (c != 0) return c;
            return Integer.compare(this.size, o.size);
        }

        @Override public boolean equals(Object obj) {
            if (this == obj) return true;
            if (!(obj instanceof EndDesc)) return false;
            EndDesc o = (EndDesc) obj;
            return this.size == o.size && this.type.compareTo(o.type) == 0;
        }

        @Override public int hashCode() {
            int h = 1;
            h = 31*h + size;
            h = 31*h + type.hashCode();
            return h;
        }
    }

    private static final class OutcomeKey implements Comparable<OutcomeKey> {
        final EndDesc e1;
        final EndDesc e2;
        final EndDesc e3; // null for pair outcomes, non-null for triplet outcomes

        final SignatureUtils.VertexSignature[] originalTypesSortedMultiset;
        final int[] originalClassSizesSortedMultiset;
        final long[] quantizedSpectrumSorted; // FIX: long[]

        OutcomeKey(EndDesc a, EndDesc b,
                   SignatureUtils.VertexSignature[] originalTypesSortedMultiset,
                   int[] originalClassSizesSortedMultiset,
                   long[] quantizedSpectrumSorted)
        {
            // canonicalize endpoints
            if (a.compareTo(b) <= 0) { this.e1 = a; this.e2 = b; }
            else { this.e1 = b; this.e2 = a; }
            this.e3 = null;

            this.originalTypesSortedMultiset = originalTypesSortedMultiset;
            this.originalClassSizesSortedMultiset = originalClassSizesSortedMultiset;
            this.quantizedSpectrumSorted = quantizedSpectrumSorted;
        }

        OutcomeKey(EndDesc a, EndDesc b, EndDesc c,
                   SignatureUtils.VertexSignature[] originalTypesSortedMultiset,
                   int[] originalClassSizesSortedMultiset,
                   long[] quantizedSpectrumSorted)
        {
            // canonicalize endpoints (sort a,b,c)
            EndDesc[] arr = new EndDesc[]{a,b,c};
            Arrays.sort(arr);
            this.e1 = arr[0];
            this.e2 = arr[1];
            this.e3 = arr[2];

            this.originalTypesSortedMultiset = originalTypesSortedMultiset;
            this.originalClassSizesSortedMultiset = originalClassSizesSortedMultiset;
            this.quantizedSpectrumSorted = quantizedSpectrumSorted;
        }

        private boolean isTriplet() { return e3 != null; }

        @Override
        public int compareTo(OutcomeKey o) {
            // Compare endpoints
            int c = this.e1.compareTo(o.e1);
            if (c != 0) return c;
            c = this.e2.compareTo(o.e2);
            if (c != 0) return c;

            // Pair keys sort before triplet keys if endpoints match
            if (this.isTriplet() != o.isTriplet()) return this.isTriplet() ? 1 : -1;

            if (this.isTriplet()) {
                c = this.e3.compareTo(o.e3);
                if (c != 0) return c;
            }

            // Compare spectrum summary (FIX: long[] compare)
            c = compareLongArrays(this.quantizedSpectrumSorted, o.quantizedSpectrumSorted);
            if (c != 0) return c;

            // Compare original-only partition summaries
            c = compareSigArrays(this.originalTypesSortedMultiset, o.originalTypesSortedMultiset);
            if (c != 0) return c;

            return compareIntArrays(this.originalClassSizesSortedMultiset, o.originalClassSizesSortedMultiset);
        }

        @Override
        public boolean equals(Object obj) {
            if (this == obj) return true;
            if (!(obj instanceof OutcomeKey)) return false;
            OutcomeKey o = (OutcomeKey) obj;

            if (!this.e1.equals(o.e1) || !this.e2.equals(o.e2)) return false;
            if (this.isTriplet() != o.isTriplet()) return false;
            if (this.isTriplet() && !this.e3.equals(o.e3)) return false;

            if (!Arrays.equals(this.quantizedSpectrumSorted, o.quantizedSpectrumSorted)) return false;
            if (!sigArrayEquals(this.originalTypesSortedMultiset, o.originalTypesSortedMultiset)) return false;
            return Arrays.equals(this.originalClassSizesSortedMultiset, o.originalClassSizesSortedMultiset);
        }

        @Override
        public int hashCode() {
            int h = 1;
            h = 31*h + e1.hashCode();
            h = 31*h + e2.hashCode();
            h = 31*h + (e3 == null ? 0 : e3.hashCode());
            h = 31*h + Arrays.hashCode(quantizedSpectrumSorted); // long[] overload
            h = 31*h + sigArrayHash(originalTypesSortedMultiset);
            h = 31*h + Arrays.hashCode(originalClassSizesSortedMultiset);
            return h;
        }
    }

    private static int compareIntArrays(int[] a, int[] b) {
        int n = Math.min(a.length, b.length);
        for (int i = 0; i < n; i++) {
            int c = Integer.compare(a[i], b[i]);
            if (c != 0) return c;
        }
        return Integer.compare(a.length, b.length);
    }

    private static int compareLongArrays(long[] a, long[] b) {
        int n = Math.min(a.length, b.length);
        for (int i = 0; i < n; i++) {
            int c = Long.compare(a[i], b[i]);
            if (c != 0) return c;
        }
        return Integer.compare(a.length, b.length);
    }

    private static int compareSigArrays(SignatureUtils.VertexSignature[] a, SignatureUtils.VertexSignature[] b) {
        int n = Math.min(a.length, b.length);
        for (int i = 0; i < n; i++) {
            int c = a[i].compareTo(b[i]);
            if (c != 0) return c;
        }
        return Integer.compare(a.length, b.length);
    }

    private static boolean sigArrayEquals(SignatureUtils.VertexSignature[] a, SignatureUtils.VertexSignature[] b) {
        if (a.length != b.length) return false;
        for (int i = 0; i < a.length; i++) {
            if (a[i].compareTo(b[i]) != 0) return false;
        }
        return true;
    }

    private static int sigArrayHash(SignatureUtils.VertexSignature[] a) {
        int h = 1;
        for (SignatureUtils.VertexSignature x : a) h = 31*h + x.hashCode();
        return h;
    }

    private static final class PairOutcomeKey implements Comparable<PairOutcomeKey> {
        final OutcomeKey first;
        final OutcomeKey second;

        PairOutcomeKey(OutcomeKey a, OutcomeKey b) {
            // IMPORTANT: canonicalize ordering so {A,B} matches {B,A}
            if (a.compareTo(b) <= 0) { this.first = a; this.second = b; }
            else { this.first = b; this.second = a; }
        }

        @Override public int compareTo(PairOutcomeKey o) {
            int c = this.first.compareTo(o.first);
            if (c != 0) return c;
            return this.second.compareTo(o.second);
        }

        @Override public boolean equals(Object obj) {
            if (this == obj) return true;
            if (!(obj instanceof PairOutcomeKey)) return false;
            PairOutcomeKey o = (PairOutcomeKey) obj;
            return this.first.equals(o.first) && this.second.equals(o.second);
        }

        @Override public int hashCode() {
            int h = 1;
            h = 31*h + first.hashCode();
            h = 31*h + second.hashCode();
            return h;
        }
    }

    /** Comparable wrapper for the multiset of probe outcomes. */
    private static final class ProbeProfile implements Comparable<ProbeProfile> {
        final NavigableMap<PairOutcomeKey, Integer> outcomeCounts;

        ProbeProfile(NavigableMap<PairOutcomeKey, Integer> outcomeCounts) {
            this.outcomeCounts = outcomeCounts;
        }

        @Override
        public int compareTo(ProbeProfile o) {
            var itA = this.outcomeCounts.entrySet().iterator();
            var itB = o.outcomeCounts.entrySet().iterator();

            while (itA.hasNext() && itB.hasNext()) {
                var ea = itA.next();
                var eb = itB.next();

                int c = ea.getKey().compareTo(eb.getKey());
                if (c != 0) return c;

                c = Integer.compare(ea.getValue(), eb.getValue());
                if (c != 0) return c;
            }

            if (itA.hasNext()) return 1;
            if (itB.hasNext()) return -1;
            return 0;
        }

        @Override
        public boolean equals(Object obj) {
            if (this == obj) return true;
            if (!(obj instanceof ProbeProfile)) return false;
            return this.outcomeCounts.equals(((ProbeProfile) obj).outcomeCounts);
        }

        @Override
        public int hashCode() {
            return outcomeCounts.hashCode();
        }
    }

    private static long choose2(long n) {
        if (n < 2) return 0L;
        return (n * (n - 1)) / 2L;
    }

    private static void probeProgressGlobal(String tag, long done, long total) {
        if (!DEBUG) return;
        if (total <= 0) return;

        long now = System.currentTimeMillis();
        long last = LAST_PROBE_PROGRESS_LOG_MS.get();
        if (now - last >= PROBE_PROGRESS_LOG_EVERY_MS) {
            if (LAST_PROBE_PROGRESS_LOG_MS.compareAndSet(last, now)) {
                long remaining = Math.max(0L, total - done);
                long pct = (total == 0) ? 100 : (done * 100L) / total;
                System.out.println(tag + " done=" + done + "/" + total + " remaining=" + remaining + " (" + pct + "%)");
            }
        }
    }

    private static ProbeProfile buildProbeProfilePairsStrongTwoOrientations(
            Graph<Integer, DefaultEdge> base,
            int u,
            int K,
            int pathLenU,
            int pathLenV,
            int nOriginal,
            AtomicLong globalDone,
            long globalTotal,
            String tag
    ) {
        TreeMap<PairOutcomeKey, Integer> outcomeCounts = new TreeMap<>();

        final int uCliqueA = K + 1;
        final int vCliqueA = K + 2;

        final int uCliqueB = K + 2; // swapped
        final int vCliqueB = K + 1;

        for (int v = 0; v < nOriginal; v++) {
            if (v == u) continue;

            OutcomeKey keyA = computeOutcomeKeyPair(
                    base, u, v, K, uCliqueA, vCliqueA, pathLenU, pathLenV, nOriginal
            );

            OutcomeKey keyB = computeOutcomeKeyPair(
                    base, u, v, K, uCliqueB, vCliqueB, pathLenV, pathLenU, nOriginal
            );

            PairOutcomeKey pairKey = new PairOutcomeKey(keyA, keyB);
            outcomeCounts.merge(pairKey, 1, Integer::sum);

            long done = globalDone.incrementAndGet();     // 1 tick per (u,v)
            probeProgressGlobal(tag, done, globalTotal);  // time-gated globally
        }

        return new ProbeProfile(outcomeCounts);
    }

    private static OutcomeKey computeOutcomeKeyPair(
            Graph<Integer, DefaultEdge> base,
            int u,
            int v,
            int K,
            int uPrivateClique,
            int vPrivateClique,
            int uPathLen,
            int vPathLen,
            int nOriginal
    ) {
        Graph<Integer, DefaultEdge> baseCopy = GraphUtils.copyGraphSameType(base);

        Graph<Integer, DefaultEdge> h = buildStrongProbeGraphPair(baseCopy, u, v, K, uPrivateClique, vPrivateClique, uPathLen, vPathLen);
        PartitionInfo info = computePartitionInfoOriginalOnly(h, nOriginal);

        SigAndSize uInfo = findSignatureAndSizeOfVertex(info.sigToVertsOriginalOnly, u);
        SigAndSize vInfo = findSignatureAndSizeOfVertex(info.sigToVertsOriginalOnly, v);

        EndDesc eu = new EndDesc(uInfo.sig, uInfo.size);
        EndDesc ev = new EndDesc(vInfo.sig, vInfo.size);

        return new OutcomeKey(
                eu, ev,
                info.originalVertexTypesSortedMultiset,
                info.originalClassSizesSortedMultiset,
                info.quantizedSpectrumSorted
        );
    }

    private static ProbeProfile buildProbeProfileTripletsStrongTwoOrientationsVW(
            Graph<Integer, DefaultEdge> base,
            int u,
            int K,
            int pathLenU,
            int pathLenV,
            int nOriginal,
            AtomicLong globalDone,
            long globalTotal,
            String tag
    ) {
        TreeMap<PairOutcomeKey, Integer> outcomeCounts = new TreeMap<>();

        // Distinct gadgets:
        final int uClique = K + 1;
        final int vClique = K + 2;
        final int wClique = K + 3;

        // Ensure pathLenW is distinct from BOTH pathLenU and pathLenV (deterministically).
        final int pathLenW = chooseThirdDistinctNonNegativePathLen(pathLenU, pathLenV);

        // Iterate unordered pairs (v<w), v!=u, w!=u
        for (int v = 0; v < nOriginal; v++) {
            if (v == u) continue;
            for (int w = v + 1; w < nOriginal; w++) {
                if (w == u) continue;

                // Orientation A: (v gets vClique/pathLenV), (w gets wClique/pathLenW)
                OutcomeKey keyA = computeOutcomeKeyTriplet(
                        base, u, v, w, K,
                        uClique, vClique, wClique,
                        pathLenU, pathLenV, pathLenW,
                        nOriginal
                );

                // Orientation B: swap v and w private gadgets + their path lengths
                OutcomeKey keyB = computeOutcomeKeyTriplet(
                        base, u, v, w, K,
                        uClique, wClique, vClique,
                        pathLenU, pathLenW, pathLenV,
                        nOriginal
                );

                PairOutcomeKey pairKey = new PairOutcomeKey(keyA, keyB);
                outcomeCounts.merge(pairKey, 1, Integer::sum);

                long done = globalDone.incrementAndGet();     // 1 tick per (u,v,w) choice
                probeProgressGlobal(tag, done, globalTotal);  // time-gated globally
            }
        }

        return new ProbeProfile(outcomeCounts);
    }

    private static int chooseThirdDistinctNonNegativePathLen(int a, int b) {
        // Deterministically pick the smallest non-negative integer not equal to a or b.
        int x = 0;
        while (x == a || x == b) x++;
        return x;
    }

    private static OutcomeKey computeOutcomeKeyTriplet(
            Graph<Integer, DefaultEdge> base,
            int u,
            int v,
            int w,
            int K,
            int uPrivateClique,
            int vPrivateClique,
            int wPrivateClique,
            int uPathLen,
            int vPathLen,
            int wPathLen,
            int nOriginal
    ) {
        Graph<Integer, DefaultEdge> baseCopy = GraphUtils.copyGraphSameType(base);

        Graph<Integer, DefaultEdge> h = buildStrongProbeGraphTriplet(
                baseCopy, u, v, w, K,
                uPrivateClique, vPrivateClique, wPrivateClique,
                uPathLen, vPathLen, wPathLen
        );

        PartitionInfo info = computePartitionInfoOriginalOnly(h, nOriginal);

        SigAndSize uInfo = findSignatureAndSizeOfVertex(info.sigToVertsOriginalOnly, u);
        SigAndSize vInfo = findSignatureAndSizeOfVertex(info.sigToVertsOriginalOnly, v);
        SigAndSize wInfo = findSignatureAndSizeOfVertex(info.sigToVertsOriginalOnly, w);

        EndDesc eu = new EndDesc(uInfo.sig, uInfo.size);
        EndDesc ev = new EndDesc(vInfo.sig, vInfo.size);
        EndDesc ew = new EndDesc(wInfo.sig, wInfo.size);

        return new OutcomeKey(
                eu, ev, ew,
                info.originalVertexTypesSortedMultiset,
                info.originalClassSizesSortedMultiset,
                info.quantizedSpectrumSorted
        );
    }

    public static List<List<Integer>> partitionVerticesBySharedCliqueProbing(
            Graph<Integer, DefaultEdge> g,
            List<Integer> verts,
            int nOriginal,
            int K,
            int pathLenU,
            int pathLenV,
            boolean debug,
            int topOutcomesToShow
    ) {
        GraphUtils.requireCanonicalZeroToNMinus1(g);
        if (K <= 0) throw new IllegalArgumentException("K must be > 0");
        if (verts == null) throw new IllegalArgumentException("verts must not be null");
        if (topOutcomesToShow < 0) throw new IllegalArgumentException("topOutcomesToShow must be >= 0");
        if (pathLenU < 0 || pathLenV < 0) throw new IllegalArgumentException("path lengths must be >= 0");
        if (pathLenU == pathLenV) throw new IllegalArgumentException("pathLenU and pathLenV must be different");

        for (int u : verts) {
            if (!g.containsVertex(u)) throw new IllegalArgumentException("vertex not in graph: " + u);
            if (u < 0 || u >= nOriginal) throw new IllegalArgumentException("vertex not in ORIGINAL range [0,nOriginal): " + u);
        }

        TreeMap<ProbeProfile, List<Integer>> profileToVerts = new TreeMap<>();
        TreeMap<Integer, ProbeProfile> perVertexProfile = new TreeMap<>();

        ExecutorService pool = Executors.newFixedThreadPool(N_THREADS);
        ConcurrentHashMap<Integer, ProbeProfile> tmpProfiles = new ConcurrentHashMap<>();
        List<Future<?>> futures = new ArrayList<>();

        // GLOBAL progress: total pair probes = |verts| * (nOriginal - 1)
        final AtomicLong globalDone = new AtomicLong(0);
        final long globalTotal = (long) verts.size() * (long) Math.max(0, nOriginal - 1);
        progressNow("[PAIR] global total=" + globalTotal + " (single-graph)");

        for (int u : verts) {
            futures.add(pool.submit(() -> {
                ProbeProfile profile = buildProbeProfilePairsStrongTwoOrientations(
                        g, u, K, pathLenU, pathLenV, nOriginal,
                        globalDone, globalTotal, "[PAIR]"
                );
                tmpProfiles.put(u, profile);
            }));
        }

        for (Future<?> f : futures) {
            try { f.get(); } catch (Exception e) { throw new RuntimeException(e); }
        }
        pool.shutdown();

        progressNow("[PAIR] global done=" + globalDone.get() + "/" + globalTotal + " (single-graph)");

        for (int u : verts) {
            ProbeProfile profile = tmpProfiles.get(u);
            perVertexProfile.put(u, profile);
            profileToVerts.computeIfAbsent(profile, __ -> new ArrayList<>()).add(u);
        }

        List<List<Integer>> groups = new ArrayList<>();
        for (var e : profileToVerts.entrySet()) groups.add(e.getValue());

        if (debug) {
            debugPrintPartition(g, verts, K, pathLenU, pathLenV, nOriginal, groups, perVertexProfile, topOutcomesToShow);
        }

        return groups;
    }

    public static TwoGraphPartitionResult partitionVerticesBySharedCliqueProbingTwoGraphs(
            Graph<Integer, DefaultEdge> g1,
            List<Integer> verts1,
            int nOriginal1,
            Graph<Integer, DefaultEdge> g2,
            List<Integer> verts2,
            int nOriginal2,
            int K,
            int pathLenU,
            int pathLenV,
            boolean debug,
            int topOutcomesToShow
    ) {
        GraphUtils.requireCanonicalZeroToNMinus1(g1);
        GraphUtils.requireCanonicalZeroToNMinus1(g2);
        if (K <= 0) throw new IllegalArgumentException("K must be > 0");
        if (verts1 == null) throw new IllegalArgumentException("verts1 must not be null");
        if (verts2 == null) throw new IllegalArgumentException("verts2 must not be null");
        if (topOutcomesToShow < 0) throw new IllegalArgumentException("topOutcomesToShow must be >= 0");
        if (pathLenU < 0 || pathLenV < 0) throw new IllegalArgumentException("path lengths must be >= 0");
        if (pathLenU == pathLenV) throw new IllegalArgumentException("pathLenU and pathLenV must be different");

        for (int u : verts1) {
            if (!g1.containsVertex(u)) throw new IllegalArgumentException("vertex not in g1: " + u);
            if (u < 0 || u >= nOriginal1) throw new IllegalArgumentException("vertex not in ORIGINAL range [0,nOriginal1): " + u);
        }
        for (int u : verts2) {
            if (!g2.containsVertex(u)) throw new IllegalArgumentException("vertex not in g2: " + u);
            if (u < 0 || u >= nOriginal2) throw new IllegalArgumentException("vertex not in ORIGINAL range [0,nOriginal2): " + u);
        }

        // First: PAIR probing
        TwoGraphPartitionResult pairs = partitionTwoGraphsInternal(
                g1, verts1, nOriginal1,
                g2, verts2, nOriginal2,
                K, pathLenU, pathLenV,
                debug, topOutcomesToShow,
                /*useTriplets*/ false
        );

        // If pair probing failed to refine (single group on BOTH sides), optionally fallback to triplets.
        boolean pairFailedToRefine =
                (pairs.groupsG1.size() == 1 && pairs.groupsG2.size() == 1
                        && verts1.size() >= 2 && verts2.size() >= 2);

        boolean skipTripletsThisRound =
                SKIP_TRIPLETS_ON_FIRST_ROUND && (CURRENT_ROUND == 1);

        if (USE_TRIPLETS_ON_PAIR_FAILURE && pairFailedToRefine && !skipTripletsThisRound) {
            if (debug) {
                System.out.println("==== Pair probing produced 1 group; falling back to TRIPLET probing ====");
                System.out.println();
            }
            return partitionTwoGraphsInternal(
                    g1, verts1, nOriginal1,
                    g2, verts2, nOriginal2,
                    K, pathLenU, pathLenV,
                    debug, topOutcomesToShow,
                    /*useTriplets*/ true
            );
        }

        if (USE_TRIPLETS_ON_PAIR_FAILURE && pairFailedToRefine && skipTripletsThisRound && debug) {
            System.out.println("==== Pair probing produced 1 group; SKIPPING triplets because CURRENT_ROUND == 1 ====");
            System.out.println();
        }

        return pairs;
    }

    private static TwoGraphPartitionResult partitionTwoGraphsInternal(
            Graph<Integer, DefaultEdge> g1,
            List<Integer> verts1,
            int nOriginal1,
            Graph<Integer, DefaultEdge> g2,
            List<Integer> verts2,
            int nOriginal2,
            int K,
            int pathLenU,
            int pathLenV,
            boolean debug,
            int topOutcomesToShow,
            boolean useTriplets
    ) {
        TreeMap<ProbeProfile, List<Integer>> profileToVerts1 = new TreeMap<>();
        TreeMap<ProbeProfile, List<Integer>> profileToVerts2 = new TreeMap<>();

        TreeMap<Integer, ProbeProfile> perVertexProfile1 = new TreeMap<>();
        TreeMap<Integer, ProbeProfile> perVertexProfile2 = new TreeMap<>();

        ExecutorService pool = Executors.newFixedThreadPool(N_THREADS);
        ConcurrentHashMap<Integer, ProbeProfile> tmp1 = new ConcurrentHashMap<>();
        ConcurrentHashMap<Integer, ProbeProfile> tmp2 = new ConcurrentHashMap<>();
        List<Future<?>> futures = new ArrayList<>();

        // GLOBAL progress across BOTH graphs
        final AtomicLong globalDone = new AtomicLong(0);

        final long total1 = useTriplets
                ? (long) verts1.size() * choose2((long) nOriginal1 - 1L)
                : (long) verts1.size() * (long) Math.max(0, nOriginal1 - 1);

        final long total2 = useTriplets
                ? (long) verts2.size() * choose2((long) nOriginal2 - 1L)
                : (long) verts2.size() * (long) Math.max(0, nOriginal2 - 1);

        final long globalTotal = total1 + total2;

        final String tag = useTriplets ? "[TRIPLET]" : "[PAIR]";
        progressNow(tag + " global total=" + globalTotal + " (G1=" + total1 + ", G2=" + total2 + ")");

        for (int u : verts1) {
            futures.add(pool.submit(() -> {
                ProbeProfile p = useTriplets
                        ? buildProbeProfileTripletsStrongTwoOrientationsVW(
                                g1, u, K, pathLenU, pathLenV, nOriginal1,
                                globalDone, globalTotal, tag
                        )
                        : buildProbeProfilePairsStrongTwoOrientations(
                                g1, u, K, pathLenU, pathLenV, nOriginal1,
                                globalDone, globalTotal, tag
                        );
                tmp1.put(u, p);
            }));
        }

        for (int u : verts2) {
            futures.add(pool.submit(() -> {
                ProbeProfile p = useTriplets
                        ? buildProbeProfileTripletsStrongTwoOrientationsVW(
                                g2, u, K, pathLenU, pathLenV, nOriginal2,
                                globalDone, globalTotal, tag
                        )
                        : buildProbeProfilePairsStrongTwoOrientations(
                                g2, u, K, pathLenU, pathLenV, nOriginal2,
                                globalDone, globalTotal, tag
                        );
                tmp2.put(u, p);
            }));
        }

        for (Future<?> f : futures) {
            try { f.get(); } catch (Exception e) { throw new RuntimeException(e); }
        }
        pool.shutdown();

        progressNow(tag + " global done=" + globalDone.get() + "/" + globalTotal);

        for (int u : verts1) {
            ProbeProfile p = tmp1.get(u);
            perVertexProfile1.put(u, p);
            profileToVerts1.computeIfAbsent(p, __ -> new ArrayList<>()).add(u);
        }

        for (int u : verts2) {
            ProbeProfile p = tmp2.get(u);
            perVertexProfile2.put(u, p);
            profileToVerts2.computeIfAbsent(p, __ -> new ArrayList<>()).add(u);
        }

        List<List<Integer>> groups1 = new ArrayList<>();
        List<ProbeProfile> groupProfiles1 = new ArrayList<>();
        for (var e : profileToVerts1.entrySet()) {
            groups1.add(e.getValue());
            groupProfiles1.add(e.getKey());
        }

        List<List<Integer>> groups2 = new ArrayList<>();
        List<ProbeProfile> groupProfiles2 = new ArrayList<>();
        for (var e : profileToVerts2.entrySet()) {
            groups2.add(e.getValue());
            groupProfiles2.add(e.getKey());
        }

        TreeMap<ProbeProfile, Integer> profileToGroupIndex2 = new TreeMap<>();
        for (int j = 0; j < groupProfiles2.size(); j++) {
            profileToGroupIndex2.put(groupProfiles2.get(j), j);
        }

        int[] mapG1toG2 = new int[groups1.size()];
        Arrays.fill(mapG1toG2, -1);
        for (int i = 0; i < groupProfiles1.size(); i++) {
            Integer j = profileToGroupIndex2.get(groupProfiles1.get(i));
            if (j != null) mapG1toG2[i] = j;
        }

        TwoGraphPartitionResult result = new TwoGraphPartitionResult(
                groups1, groups2, mapG1toG2, groupProfiles1, groupProfiles2
        );

        if (debug) {
            System.out.println(useTriplets
                    ? "==== Two-graph probing debug (TRIPLETS; 2 orientations swapping v/w gadgets) ===="
                    : "==== Two-graph probing debug (PAIRS; 2 orientations per pair) ====");
            System.out.println("K=" + K + "  pathLenU=" + pathLenU + "  pathLenV=" + pathLenV);
            System.out.println("G1 candidates=" + verts1);
            System.out.println("G2 candidates=" + verts2);
            System.out.println();

            System.out.println("---- G1 groups (" + groups1.size() + ") ----");
            debugPrintPartition(g1, verts1, K, pathLenU, pathLenV, nOriginal1, groups1, perVertexProfile1, topOutcomesToShow);

            System.out.println("---- G2 groups (" + groups2.size() + ") ----");
            debugPrintPartition(g2, verts2, K, pathLenU, pathLenV, nOriginal2, groups2, perVertexProfile2, topOutcomesToShow);

            System.out.println("==== Cross-graph group mapping (by identical ProbeProfile) ====");
            int mapped = 0;
            for (int i = 0; i < result.groupIndexG1ToG2.length; i++) {
                int j = result.groupIndexG1ToG2[i];
                if (j >= 0) mapped++;
                System.out.println("G1 group " + i + " -> " + (j >= 0 ? ("G2 group " + j) : "UNMAPPED"));
            }
            System.out.println("Mapped " + mapped + " / " + result.groupIndexG1ToG2.length + " G1 groups");
            System.out.println();
        }

        return result;
    }

    public static final class TwoGraphPartitionResult {
        public final List<List<Integer>> groupsG1;
        public final List<List<Integer>> groupsG2;

        /** map[i] = j means groupsG1[i] matches groupsG2[j] by identical ProbeProfile; -1 means unmapped. */
        public final int[] groupIndexG1ToG2;

        public final List<ProbeProfile> groupProfilesG1;
        public final List<ProbeProfile> groupProfilesG2;

        public TwoGraphPartitionResult(
                List<List<Integer>> groupsG1,
                List<List<Integer>> groupsG2,
                int[] groupIndexG1ToG2,
                List<ProbeProfile> groupProfilesG1,
                List<ProbeProfile> groupProfilesG2
        ) {
            this.groupsG1 = groupsG1;
            this.groupsG2 = groupsG2;
            this.groupIndexG1ToG2 = groupIndexG1ToG2;
            this.groupProfilesG1 = groupProfilesG1;
            this.groupProfilesG2 = groupProfilesG2;
        }
    }

    private static void debugPrintPartition(
            Graph<Integer, DefaultEdge> g,
            List<Integer> verts,
            int K,
            int pathLenU,
            int pathLenV,
            int nOriginal,
            List<List<Integer>> groups,
            TreeMap<Integer, ProbeProfile> perVertexProfile,
            int topOutcomesToShow
    ) {
        System.out.println("==== Probing debug ====");
        System.out.println("K=" + K + "  pathLenU=" + pathLenU + "  pathLenV=" + pathLenV + "  candidates=" + verts);
        System.out.println("Groups (" + groups.size() + "): " + groups);
        System.out.println();

        for (int i = 0; i < groups.size(); i++) {
            List<Integer> group = groups.get(i);
            int rep = group.get(0);
            ProbeProfile prof = perVertexProfile.get(rep);

            System.out.println("Group " + (i + 1) + " size=" + group.size() + " members=" + group);
            System.out.println("  representative=" + rep);
            System.out.println("  distinct probe outcomes=" + prof.outcomeCounts.size());
            if (topOutcomesToShow > 0) {
                printTopOutcomes(prof.outcomeCounts, topOutcomesToShow);
            }
            System.out.println();
        }
    }

    private static void printTopOutcomes(NavigableMap<PairOutcomeKey, Integer> outcomeCounts, int top) {
        List<java.util.Map.Entry<PairOutcomeKey, Integer>> list = new ArrayList<>(outcomeCounts.entrySet());
        list.sort((a, b) -> {
            int c = Integer.compare(b.getValue(), a.getValue());
            if (c != 0) return c;
            return a.getKey().compareTo(b.getKey());
        });

        int shown = 0;
        for (var e : list) {
            if (shown++ >= top) break;
            System.out.println("    count=" + e.getValue() + "  pairOutcomeKeyHash=" + e.getKey().hashCode());
        }
        if (list.size() > top) {
            System.out.println("    ... (" + (list.size() - top) + " more distinct outcomes)");
        }
    }

    private static boolean isIsomorphismOnRawGraphs(
            Graph<Integer, DefaultEdge> g1Raw,
            Graph<Integer, DefaultEdge> g2Raw,
            int[] map,
            int nOriginal
    ) {
        for (int u = 0; u < nOriginal; u++) {
            for (int v = u + 1; v < nOriginal; v++) {
                boolean e1 = g1Raw.containsEdge(u, v);
                boolean e2 = g2Raw.containsEdge(map[u], map[v]);
                if (e1 != e2) return false;
            }
        }
        return true;
    }

    public static List<List<Integer>> smallestCommonOriginalEquivalenceClassVertsAtLeast2(
            Graph<Integer, DefaultEdge> g1,
            long[][] q1,
            Graph<Integer, DefaultEdge> g2,
            long[][] q2,
            int nOriginal,
            int bfsLayers
    ) {
        SignatureUtils.vertexEquivalenceClassesOriginalOnly(g1, q1, nOriginal, bfsLayers);
        SignatureUtils.vertexEquivalenceClassesOriginalOnly(g2, q2, nOriginal, bfsLayers);

        NavigableMap<SignatureUtils.VertexSignature, List<Integer>> t1 =
                SignatureUtils.signatureToVerticesOriginalOnly(g1, q1, nOriginal, bfsLayers);
        NavigableMap<SignatureUtils.VertexSignature, List<Integer>> t2 =
                SignatureUtils.signatureToVerticesOriginalOnly(g2, q2, nOriginal, bfsLayers);

        SignatureUtils.VertexSignature bestKey = null;
        int bestSize = Integer.MAX_VALUE;

        for (var e1 : t1.entrySet()) {
            SignatureUtils.VertexSignature key = e1.getKey();
            List<Integer> a = e1.getValue();
            if (a.size() < 2) continue;

            List<Integer> b = t2.get(key);
            if (b == null) continue;
            if (b.size() < 2) continue;
            if (a.size() != b.size()) continue;

            int sz = a.size();
            if (sz < bestSize || (sz == bestSize && (bestKey == null || key.compareTo(bestKey) < 0))) {
                bestSize = sz;
                bestKey = key;
            }
        }

        if (bestKey == null) return null;

        List<Integer> out1 = new ArrayList<>(t1.get(bestKey));
        List<Integer> out2 = new ArrayList<>(t2.get(bestKey));
        out1.sort(Integer::compareTo);
        out2.sort(Integer::compareTo);

        return Arrays.asList(out1, out2);
    }

    private static void appendCsvLineQuiet(Path logPath, String graphPath, String status, String message) {
        try {
            appendCsvLine(logPath, graphPath, status, message);
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
    }

    private static void appendCsvLine(Path logPath, String graphPath, String status, String message) throws IOException {
        // Open in append mode; no explicit file locking.
        try (BufferedWriter w = Files.newBufferedWriter(
                logPath,
                StandardCharsets.UTF_8,
                StandardOpenOption.CREATE,
                StandardOpenOption.WRITE,
                StandardOpenOption.APPEND
        )) {
            w.write(csvEscape(graphPath));
            w.write(',');
            w.write(csvEscape(status));
            w.write(',');
            w.write(csvEscape(message));
            w.newLine();
            w.flush();
        }
    }

    private static Path deriveLogPath(Path baseLog, String suffix) {
        String name = baseLog.getFileName().toString();

        int dot = name.lastIndexOf('.');
        String stem = (dot >= 0) ? name.substring(0, dot) : name;
        String ext  = (dot >= 0) ? name.substring(dot) : ".csv";

        return baseLog.resolveSibling(stem + "_" + suffix + ext);
    }

    private static void ensureCsvHeader(Path logPath) throws IOException {
        boolean needsHeader = (!Files.exists(logPath)) || (Files.size(logPath) == 0);
        if (needsHeader) {
            appendCsvLine(logPath, "graph_path", "status", "message");
        }
    }

    private static String csvEscape(String s) {
        if (s == null) return "";
        boolean needQuotes =
                s.indexOf(',') >= 0 || s.indexOf('"') >= 0 || s.indexOf('\n') >= 0 || s.indexOf('\r') >= 0;
        if (!needQuotes) return s;
        return "\"" + s.replace("\"", "\"\"") + "\"";
    }

    private static String safeMessage(Throwable t) {
        String msg = t.getMessage();
        if (msg != null && !msg.isBlank()) return msg;
        return t.getClass().getName();
    }

    private static TwinKind twinKindInRaw(Graph<Integer, DefaultEdge> raw, List<Integer> verts, int nOriginal) {
        if (verts == null || verts.size() <= 1) return TwinKind.NONE;

        int r = verts.get(0);

        // Determine kind by adjacency between r and the second element.
        int v0 = verts.get(1);
        boolean rAdjV0 = raw.containsEdge(r, v0);
        TwinKind kind = rAdjV0 ? TwinKind.TRUE_TWINS : TwinKind.FALSE_TWINS;

        for (int idx = 1; idx < verts.size(); idx++) {
            int v = verts.get(idx);

            // Must match adjacency kind with representative
            boolean rAdjV = raw.containsEdge(r, v);
            if (kind == TwinKind.TRUE_TWINS) {
                if (!rAdjV) return TwinKind.NONE;
            } else { // FALSE_TWINS
                if (rAdjV) return TwinKind.NONE;
            }

            // Neighborhood equality (excluding r and v)
            for (int x = 0; x < nOriginal; x++) {
                if (x == r || x == v) continue;

                boolean rAdjX = raw.containsEdge(r, x);
                boolean vAdjX = raw.containsEdge(v, x);
                if (rAdjX != vAdjX) return TwinKind.NONE;
            }
        }

        return kind;
    }

    private static boolean allClassesSingletonOrTwinsAndMatchAcrossGraphs(
            Graph<Integer, DefaultEdge> g1Aug, long[][] q1,
            Graph<Integer, DefaultEdge> g2Aug, long[][] q2,
            Graph<Integer, DefaultEdge> g1raw,
            Graph<Integer, DefaultEdge> g2raw,
            int nOriginal,
            int bfsLayers
    ) {
        NavigableMap<SignatureUtils.VertexSignature, List<Integer>> t1 =
                SignatureUtils.signatureToVerticesOriginalOnly(g1Aug, q1, nOriginal, bfsLayers);
        NavigableMap<SignatureUtils.VertexSignature, List<Integer>> t2 =
                SignatureUtils.signatureToVerticesOriginalOnly(g2Aug, q2, nOriginal, bfsLayers);

        if (t1 == null || t2 == null) return false;

        // Strongest check for "same classes": same key set (uses VertexSignature.compareTo/equals)
        if (t1.size() != t2.size()) return false;
        if (!t1.navigableKeySet().equals(t2.navigableKeySet())) return false;

        // Sanity: classes must partition EXACTLY the original vertex set on each side
        boolean[] seen1 = new boolean[nOriginal];
        boolean[] seen2 = new boolean[nOriginal];
        int count1 = 0;
        int count2 = 0;

        for (List<Integer> cls : t1.values()) {
            if (cls == null) return false;
            for (int v : cls) {
                if (v < 0 || v >= nOriginal) return false;
                if (seen1[v]) return false;
                seen1[v] = true;
                count1++;
            }
        }

        for (List<Integer> cls : t2.values()) {
            if (cls == null) return false;
            for (int v : cls) {
                if (v < 0 || v >= nOriginal) return false;
                if (seen2[v]) return false;
                seen2[v] = true;
                count2++;
            }
        }

        if (count1 != nOriginal || count2 != nOriginal) return false;

        // Now verify: each corresponding class is singleton OR a twin class in RAW,
        // and twin kind matches across graphs for that signature key.
        for (var e : t1.entrySet()) {
            SignatureUtils.VertexSignature key = e.getKey();
            List<Integer> a = e.getValue();
            List<Integer> b = t2.get(key);
            if (b == null) return false;

            if (a.size() != b.size()) return false;

            if (a.size() <= 1) continue; // singleton OK

            TwinKind k1 = twinKindInRaw(g1raw, a, nOriginal);
            if (k1 == TwinKind.NONE) return false;

            TwinKind k2 = twinKindInRaw(g2raw, b, nOriginal);
            if (k2 == TwinKind.NONE) return false;

            if (k1 != k2) return false;
        }

        return true;
    }

    private static int[] buildOriginalMappingByTypeAllowingTwins(
            Graph<Integer, DefaultEdge> g1Aug, long[][] q1,
            Graph<Integer, DefaultEdge> g2Aug, long[][] q2,
            Graph<Integer, DefaultEdge> g1raw,
            Graph<Integer, DefaultEdge> g2raw,
            int nOriginal,
            int bfsLayers
    ) {
        NavigableMap<SignatureUtils.VertexSignature, List<Integer>> t1 =
                SignatureUtils.signatureToVerticesOriginalOnly(g1Aug, q1, nOriginal, bfsLayers);
        NavigableMap<SignatureUtils.VertexSignature, List<Integer>> t2 =
                SignatureUtils.signatureToVerticesOriginalOnly(g2Aug, q2, nOriginal, bfsLayers);

        int[] map = new int[nOriginal];
        Arrays.fill(map, -1);

        for (var e : t1.entrySet()) {
            SignatureUtils.VertexSignature key = e.getKey();
            List<Integer> a = new ArrayList<>(e.getValue());
            List<Integer> b = t2.get(key);
            if (b == null) throw new IllegalStateException("Missing signature in G2: " + key);

            List<Integer> bb = new ArrayList<>(b);

            if (a.size() != bb.size()) {
                throw new IllegalStateException("Class size mismatch for signature " + key + ": " + a.size() + " vs " + bb.size());
            }

            // If non-singleton, require it's actually a twin-class in both raw graphs (and same kind)
            if (a.size() > 1) {
                TwinKind k1 = twinKindInRaw(g1raw, a, nOriginal);
                TwinKind k2 = twinKindInRaw(g2raw, bb, nOriginal);
                if (k1 == TwinKind.NONE || k2 == TwinKind.NONE || k1 != k2) {
                    throw new IllegalStateException("Non-singleton class is not matching twin-class for signature " + key);
                }
            }

            // Arbitrary but deterministic bijection within the class
            a.sort(Integer::compareTo);
            bb.sort(Integer::compareTo);

            for (int i = 0; i < a.size(); i++) {
                map[a.get(i)] = bb.get(i);
            }
        }

        // Validate permutation
        boolean[] seen = new boolean[nOriginal];
        for (int u = 0; u < nOriginal; u++) {
            int v = map[u];
            if (v < 0 || v >= nOriginal) throw new IllegalStateException("Mapping incomplete at u=" + u);
            if (seen[v]) throw new IllegalStateException("Mapping not injective at v=" + v);
            seen[v] = true;
        }

        return map;
    }

    private static java.util.Set<String> readFirstColumnCsv(Path csv) throws IOException {
        java.util.Set<String> out = new java.util.HashSet<>();
        if (!Files.exists(csv)) return out;

        try (Stream<String> lines = Files.lines(csv, StandardCharsets.UTF_8)) {
            lines.forEach(line -> {
                if (line == null || line.isBlank()) return;
                String first = parseFirstCsvField(line);
                if (first != null && !first.isBlank()) out.add(first);
            });
        }
        return out;
    }

    // Minimal CSV first-field parser for your format:
    // - supports quoted fields
    // - supports "" inside quotes
    private static String parseFirstCsvField(String line) {
        int i = 0;
        int n = line.length();
        if (n == 0) return "";

        if (line.charAt(0) != '"') {
            int comma = line.indexOf(',');
            return (comma >= 0) ? line.substring(0, comma) : line;
        }

        // quoted
        StringBuilder sb = new StringBuilder();
        i++; // skip opening quote
        while (i < n) {
            char c = line.charAt(i);
            if (c == '"') {
                // "" => literal "
                if (i + 1 < n && line.charAt(i + 1) == '"') {
                    sb.append('"');
                    i += 2;
                    continue;
                }
                // end quote
                i++;
                break;
            } else {
                sb.append(c);
                i++;
            }
        }
        return sb.toString();
    }

    private static boolean isTooHugeToLoad(Path graphFile) throws IOException {
        // 1) raw file byte size guard (fastest)
        long bytes = Files.size(graphFile);
        if (MAX_INPUT_BYTES > 0 && bytes > MAX_INPUT_BYTES) return true;

        // 2) quick DIMACS header scan
        DimacsHeader h = tryReadDimacsHeader(graphFile, 500);
        if (h != null) {
            if (h.n > MAX_ORIG_NODES) return true;
            if (MAX_ORIG_EDGES > 0 && h.m > MAX_ORIG_EDGES) return true;
        }

        return false;
    }

    private static final class DimacsHeader {
        final int n;
        final int m;
        DimacsHeader(int n, int m) { this.n = n; this.m = m; }
    }

    private static DimacsHeader tryReadDimacsHeader(Path graphFile, int maxLines) throws IOException {
        try (java.io.BufferedReader br = Files.newBufferedReader(graphFile, StandardCharsets.UTF_8)) {
            String line;
            int read = 0;
            while ((line = br.readLine()) != null && read++ < maxLines) {
                line = line.trim();
                if (line.isEmpty()) continue;
                if (line.startsWith("c")) continue;
                if (line.startsWith("p")) {
                    String[] parts = line.split("\\s+");
                    if (parts.length >= 4) {
                        int n, m;
                        if (parts.length >= 5) {
                            n = parseIntSafe(parts[2]);
                            m = parseIntSafe(parts[3]);
                        } else {
                            n = parseIntSafe(parts[1]);
                            m = parseIntSafe(parts[2]);
                        }
                        if (n > 0 && m >= 0) return new DimacsHeader(n, m);
                    }
                }
                if (line.startsWith("e") || line.startsWith("a")) break;
            }
        } catch (Exception ignored) {
        }
        return null;
    }

    private static int parseIntSafe(String s) {
        try { return Integer.parseInt(s); } catch (Exception e) { return -1; }
    }

    private static boolean spectralDimensionSucceeds(Graph<Integer, DefaultEdge> g) {
        try {
            int n = g.vertexSet().size();
            DMatrixRMaj V = new DMatrixRMaj(n, n);
            double[] evals = EigenUtils.computeSpectrumFromGraph(
                    g, EigenUtils.GraphMatrixType.LAPLACIAN, V
            );
            SpectralUtils.DsResult ds = SpectralUtils.estimateDs(
                    evals, Constants.SPECTRAL_EXPONENT, Constants.DIFFUSION_CONSTANT
            );
            double d = ds.ds;
            return Double.isFinite(d) && d > 0.0;
        } catch (Throwable t) {
            return false;
        }
    }

    public static void progress(String msg) {
        if (!DEBUG) return;
        long now = System.currentTimeMillis();
        long last = LAST_PROGRESS_LOG_MS.get();
        if (now - last >= PROGRESS_LOG_EVERY_MS) {
            if (LAST_PROGRESS_LOG_MS.compareAndSet(last, now)) {
                System.out.println(msg);
            }
        }
    }

    public static void progressNow(String msg) {
        if (!DEBUG) return;
        System.out.println(msg);
    }
    
 // --- CLI override support ----------------------------------------------------

    private static void applyCliOverrides(String[] args) {
        java.util.Map<String, String> opts = parseArgs(args);

        // ints
        COEFFICIENTS = optInt(opts, "COEFFICIENTS", COEFFICIENTS);
        N_THREADS = optInt(opts, "N_THREADS", N_THREADS);
        ORIG_GADGET_SIZE = optInt(opts, "ORIG_GADGET_SIZE", ORIG_GADGET_SIZE);
        ORIG_SUBD = optInt(opts, "ORIG_SUBD", ORIG_SUBD);
        MAX_ORIG_NODES = optInt(opts, "MAX_ORIG_NODES", MAX_ORIG_NODES);
        MIN_VERTICES_AFTER_SUBDIVISION = optInt(opts, "MIN_VERTICES_AFTER_SUBDIVISION", MIN_VERTICES_AFTER_SUBDIVISION);
        MAX_VERTICES_AFTER_SUBDIVISION = optInt(opts, "MAX_VERTICES_AFTER_SUBDIVISION", MAX_VERTICES_AFTER_SUBDIVISION);
        TOP_OUTCOMES_TO_SHOW = optInt(opts, "TOP_OUTCOMES_TO_SHOW", TOP_OUTCOMES_TO_SHOW);
        MAX_ORIG_EDGES = optInt(opts, "MAX_ORIG_EDGES", MAX_ORIG_EDGES);

        // longs
        SEED = optLong(opts, "SEED", SEED);
        PROGRESS_LOG_EVERY_MS = optLong(opts, "PROGRESS_LOG_EVERY_MS", PROGRESS_LOG_EVERY_MS);
        PROBE_PROGRESS_LOG_EVERY_MS = optLong(opts, "PROBE_PROGRESS_LOG_EVERY_MS", PROBE_PROGRESS_LOG_EVERY_MS);
        MAX_INPUT_BYTES = optLong(opts, "MAX_INPUT_BYTES", MAX_INPUT_BYTES);

        // booleans
        USE_TRIPLETS_ON_PAIR_FAILURE = optBool(opts, "USE_TRIPLETS_ON_PAIR_FAILURE", USE_TRIPLETS_ON_PAIR_FAILURE);
        SKIP_TRIPLETS_ON_FIRST_ROUND = optBool(opts, "SKIP_TRIPLETS_ON_FIRST_ROUND", SKIP_TRIPLETS_ON_FIRST_ROUND);
        DEBUG = optBool(opts, "DEBUG", DEBUG);

        // strings / paths
        GRAPH_PATH = optString(opts, "GRAPH_PATH", GRAPH_PATH);
        GRAPH_FOLDER = optString(opts, "GRAPH_FOLDER", GRAPH_FOLDER);
        OUTPUT_LOG = optString(opts, "OUTPUT_LOG", OUTPUT_LOG);

        // doubles
        EIGEN_DEFAULT_SCALE_OV      = optDouble(opts, "EIGEN_DEFAULT_SCALE_OV", EIGEN_DEFAULT_SCALE_OV);
        DS_EPS_OV                   = optDouble(opts, "DS_EPS_OV", DS_EPS_OV);
        CURV_COEFF_DEFAULT_SCALE_OV = optDouble(opts, "CURV_COEFF_DEFAULT_SCALE_OV", CURV_COEFF_DEFAULT_SCALE_OV);
        
        if (DEBUG && !opts.isEmpty()) {
            System.out.println("Applied CLI overrides: " + opts);
        }
    }

    private static java.util.Map<String, String> parseArgs(String[] args) {
        java.util.Map<String, String> out = new java.util.HashMap<>();
        if (args == null) return out;

        for (int i = 0; i < args.length; i++) {
            String a = args[i];
            if (a == null) continue;
            a = a.trim();
            if (!a.startsWith("--")) continue;

            // --no-DEBUG
            if (a.startsWith("--no-")) {
                String key = a.substring("--no-".length()).trim();
                if (!key.isEmpty()) out.put(key, "false");
                continue;
            }

            // --KEY=value
            int eq = a.indexOf('=');
            if (eq >= 0) {
                String key = a.substring(2, eq).trim();
                String val = a.substring(eq + 1).trim();
                if (!key.isEmpty()) out.put(key, stripQuotes(val));
                continue;
            }

            // --KEY [value]  OR bare flag: --DEBUG
            String key = a.substring(2).trim();
            if (key.isEmpty()) continue;

            // If next token is another flag or missing => treat as boolean true
            String val = "true";
            if (i + 1 < args.length) {
                String nxt = args[i + 1];
                if (nxt != null && !nxt.trim().startsWith("--")) {
                    val = nxt.trim();
                    i++;
                }
            }
            out.put(key, stripQuotes(val));
        }

        return out;
    }

    private static String stripQuotes(String s) {
        if (s == null) return null;
        s = s.trim();
        if (s.length() >= 2) {
            if ((s.startsWith("\"") && s.endsWith("\"")) || (s.startsWith("'") && s.endsWith("'"))) {
                return s.substring(1, s.length() - 1);
            }
        }
        return s;
    }

    private static int optInt(java.util.Map<String, String> m, String k, int def) {
        String v = m.get(k);
        if (v == null) return def;
        try { return Integer.parseInt(v); } catch (Exception e) { return def; }
    }

    private static double optDouble(java.util.Map<String, String> m, String k, double def) {
        String v = m.get(k);
        if (v == null) return def;
        try { return Double.parseDouble(v); } catch (Exception e) { return def; }
    }
    
    private static long optLong(java.util.Map<String, String> m, String k, long def) {
        String v = m.get(k);
        if (v == null) return def;
        try { return Long.parseLong(v); } catch (Exception e) { return def; }
    }

    private static boolean optBool(java.util.Map<String, String> m, String k, boolean def) {
        String v = m.get(k);
        if (v == null) return def;
        if (v.isEmpty()) return true;
        return v.equalsIgnoreCase("true") || v.equalsIgnoreCase("1") || v.equalsIgnoreCase("yes") || v.equalsIgnoreCase("y");
    }

    private static String optString(java.util.Map<String, String> m, String k, String def) {
        String v = m.get(k);
        return (v == null) ? def : v;
    }
}