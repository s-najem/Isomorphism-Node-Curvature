package iso;

import org.jgrapht.Graph;
import org.jgrapht.Graphs;
import org.jgrapht.graph.DefaultEdge;

import java.util.*;

public final class SignatureUtils {

    private SignatureUtils() {}

    public static String summarizeClassSizes(EquivalenceClasses classes) {
        Objects.requireNonNull(classes, "classes is null");
        int[] sizes = classes.sizes;

        // size -> count
        TreeMap<Integer, Integer> freq = new TreeMap<>();
        for (int sz : sizes) {
            freq.merge(sz, 1, Integer::sum);
        }

        StringBuilder sb = new StringBuilder();
        sb.append("nClasses=").append(classes.nClasses).append(", countsBySize={");
        boolean first = true;
        for (var e : freq.entrySet()) {
            if (!first) sb.append(", ");
            first = false;
            sb.append(e.getKey()).append(":").append(e.getValue());
        }
        sb.append("}");
        return sb.toString();
    }
    
    /**
     * Build vertex signatures only for vertices in [startInclusive, endExclusive).
     * Returned array length = endExclusive-startInclusive.
     * Each VertexSignature.root remains the original vertex id.
     */
    public static VertexSignature[] buildVertexSignaturesRange(
            Graph<Integer, DefaultEdge> g,
            long[][] qCoeffs,
            int startInclusive,
            int endExclusive,
            int bfsLayers
    ) {
        Objects.requireNonNull(g, "graph is null");
        Objects.requireNonNull(qCoeffs, "qCoeffs is null");
        if (bfsLayers <= 0) throw new IllegalArgumentException("bfsLayers must be >= 1");

        // Strong invariants (connectedness + canonical labels) for stability/invariance.
        EigenUtils.requireSimpleConnectedCanonical(g);

        int n = g.vertexSet().size();
        if (qCoeffs.length != n) {
            throw new IllegalArgumentException("qCoeffs length must equal n. qCoeffs.length=" + qCoeffs.length + " n=" + n);
        }

        if (startInclusive < 0 || endExclusive < startInclusive || endExclusive > n) {
            throw new IllegalArgumentException("Invalid range [" + startInclusive + "," + endExclusive + ") for n=" + n);
        }

        int m = endExclusive - startInclusive;
        VertexSignature[] out = new VertexSignature[m];
        for (int root = startInclusive; root < endExclusive; root++) {
            out[root - startInclusive] = buildVertexSignature(g, qCoeffs, root, bfsLayers);
        }
        return out;
    }

    /**
     * Convenience: build signatures for ORIGINAL vertices only: [0, nOriginal).
     * nOriginal is usually the original graph size BEFORE augmentations (universal vertex / cliques / subdivision).
     */
    public static VertexSignature[] buildVertexSignaturesOriginalOnly(
            Graph<Integer, DefaultEdge> g,
            long[][] qCoeffs,
            int nOriginal,
            int bfsLayers
    ) {
        int n = g.vertexSet().size();
        int end = Math.min(nOriginal, n);
        return buildVertexSignaturesRange(g, qCoeffs, 0, end, bfsLayers);
    }

    /**
     * Group vertices in [startInclusive, endExclusive) by VertexSignature.
     *
     * IMPORTANT: returns a TreeMap keyed by signature using compareTo (NOT hashing).
     */
    public static NavigableMap<VertexSignature, List<Integer>> signatureToVerticesRange(
            Graph<Integer, DefaultEdge> g,
            long[][] qCoeffs,
            int startInclusive,
            int endExclusive,
            int bfsLayers
    ) {
        VertexSignature[] sigs = buildVertexSignaturesRange(g, qCoeffs, startInclusive, endExclusive, bfsLayers);

        // Use TreeMap with compareTo to avoid relying on hash identity.
        TreeMap<VertexSignature, List<Integer>> map = new TreeMap<>(VertexSignature::compareTo);

        for (VertexSignature vs : sigs) {
            map.computeIfAbsent(vs, k -> new ArrayList<>()).add(vs.root);
        }

        // Make deterministic member ordering
        for (List<Integer> lst : map.values()) {
            Collections.sort(lst);
        }

        return map;
    }

    /** Convenience: signatures -> vertices for ORIGINAL vertices only [0, nOriginal). */
    public static NavigableMap<VertexSignature, List<Integer>> signatureToVerticesOriginalOnly(
            Graph<Integer, DefaultEdge> g,
            long[][] qCoeffs,
            int nOriginal,
            int bfsLayers
    ) {
        int n = g.vertexSet().size();
        int end = Math.min(nOriginal, n);
        return signatureToVerticesRange(g, qCoeffs, 0, end, bfsLayers);
    }

    /**
     * Equivalence classes for vertices in [startInclusive, endExclusive).
     *
     * Output convention:
     * - classIdPerVertex has length = endExclusive-startInclusive, indexing is local (vLocal = v - startInclusive)
     * - representative[] stores the GLOBAL vertex id
     */
    public static EquivalenceClasses vertexEquivalenceClassesRange(
            Graph<Integer, DefaultEdge> g,
            long[][] qCoeffs,
            int startInclusive,
            int endExclusive,
            int bfsLayers
    ) {
        VertexSignature[] vSigs = buildVertexSignaturesRange(g, qCoeffs, startInclusive, endExclusive, bfsLayers);

        int m = vSigs.length;

        Integer[] idx = new Integer[m];
        for (int i = 0; i < m; i++) idx[i] = i;

        Arrays.sort(idx, (a, b) -> vSigs[a].compareTo(vSigs[b]));

        int[] classIdLocal = new int[m];
        Arrays.fill(classIdLocal, -1);

        List<Integer> reps = new ArrayList<>();
        List<Integer> sizes = new ArrayList<>();

        int c = -1;
        int runStart = 0;
        while (runStart < m) {
            c++;
            int repLocal = idx[runStart];
            int repGlobal = vSigs[repLocal].root;
            reps.add(repGlobal);

            int runEnd = runStart + 1;
            while (runEnd < m && vSigs[idx[runStart]].compareTo(vSigs[idx[runEnd]]) == 0) {
                runEnd++;
            }

            int sz = runEnd - runStart;
            sizes.add(sz);

            for (int k = runStart; k < runEnd; k++) {
                classIdLocal[idx[k]] = c;
            }

            runStart = runEnd;
        }

        int[] sizesArr = new int[sizes.size()];
        int[] repsArr = new int[reps.size()];
        for (int i = 0; i < sizesArr.length; i++) sizesArr[i] = sizes.get(i);
        for (int i = 0; i < repsArr.length; i++) repsArr[i] = reps.get(i);

        return new EquivalenceClasses(classIdLocal, repsArr, sizesArr);
    }

    /** Convenience: equivalence classes for ORIGINAL vertices [0, nOriginal). */
    public static EquivalenceClasses vertexEquivalenceClassesOriginalOnly(
            Graph<Integer, DefaultEdge> g,
            long[][] qCoeffs,
            int nOriginal,
            int bfsLayers
    ) {
        int n = g.vertexSet().size();
        int end = Math.min(nOriginal, n);
        return vertexEquivalenceClassesRange(g, qCoeffs, 0, end, bfsLayers);
    }

    /** Build all vertex signatures for a graph (using Constants.SIGNATURE_BFS_LAYERS). */
    public static VertexSignature[] buildVertexSignatures(
            Graph<Integer, DefaultEdge> g,
            long[][] qCoeffs
    ) {
        return buildVertexSignatures(g, qCoeffs, Constants.SIGNATURE_BFS_LAYERS);
    }

    /** Build all vertex signatures for a graph with explicit BFS layers. */
    public static VertexSignature[] buildVertexSignatures(
            Graph<Integer, DefaultEdge> g,
            long[][] qCoeffs,
            int bfsLayers
    ) {
        Objects.requireNonNull(g, "graph is null");
        Objects.requireNonNull(qCoeffs, "qCoeffs is null");
        if (bfsLayers <= 0) throw new IllegalArgumentException("bfsLayers must be >= 1");

        // Strong invariants (connectedness + canonical labels) for stability/invariance.
        EigenUtils.requireSimpleConnectedCanonical(g);

        int n = g.vertexSet().size();
        if (qCoeffs.length != n) {
            throw new IllegalArgumentException("qCoeffs length must equal n. qCoeffs.length=" + qCoeffs.length + " n=" + n);
        }

        VertexSignature[] out = new VertexSignature[n];
        for (int root = 0; root < n; root++) {
            out[root] = buildVertexSignature(g, qCoeffs, root, bfsLayers);
        }
        return out;
    }

    /** Build a graph signature (multiset of vertex signatures), canonically sorted. */
    public static GraphSignature buildGraphSignature(
            Graph<Integer, DefaultEdge> g,
            long[][] qCoeffs
    ) {
        return buildGraphSignature(g, qCoeffs, Constants.SIGNATURE_BFS_LAYERS);
    }

    public static GraphSignature buildGraphSignature(
            Graph<Integer, DefaultEdge> g,
            long[][] qCoeffs,
            int bfsLayers
    ) {
        VertexSignature[] vSigs = buildVertexSignatures(g, qCoeffs, bfsLayers);

        // Canonical multiset: sort vertex signatures ignoring vertex ids.
        VertexSignature[] sorted = vSigs.clone();
        Arrays.sort(sorted);

        return new GraphSignature(sorted);
    }

    /** Compute vertex-signature equivalence classes for a graph. */
    public static EquivalenceClasses vertexEquivalenceClasses(
            Graph<Integer, DefaultEdge> g,
            long[][] qCoeffs
    ) {
        return vertexEquivalenceClasses(g, qCoeffs, Constants.SIGNATURE_BFS_LAYERS);
    }

    public static EquivalenceClasses vertexEquivalenceClasses(
            Graph<Integer, DefaultEdge> g,
            long[][] qCoeffs,
            int bfsLayers
    ) {
        VertexSignature[] vSigs = buildVertexSignatures(g, qCoeffs, bfsLayers);

        // Sort indices by signature to group equals
        int n = vSigs.length;
        Integer[] idx = new Integer[n];
        for (int i = 0; i < n; i++) idx[i] = i;

        Arrays.sort(idx, (a, b) -> vSigs[a].compareTo(vSigs[b]));

        int[] classId = new int[n];
        Arrays.fill(classId, -1);

        List<Integer> reps = new ArrayList<>();
        List<Integer> sizes = new ArrayList<>();

        int c = -1;
        int runStart = 0;

        while (runStart < n) {
            c++;
            int repVertex = idx[runStart];
            reps.add(repVertex);

            int runEnd = runStart + 1;
            while (runEnd < n && vSigs[idx[runStart]].compareTo(vSigs[idx[runEnd]]) == 0) {
                runEnd++;
            }

            int sz = runEnd - runStart;
            sizes.add(sz);

            for (int k = runStart; k < runEnd; k++) {
                classId[idx[k]] = c;
            }

            runStart = runEnd;
        }

        int[] sizesArr = new int[sizes.size()];
        int[] repsArr = new int[reps.size()];
        for (int i = 0; i < sizesArr.length; i++) sizesArr[i] = sizes.get(i);
        for (int i = 0; i < repsArr.length; i++) repsArr[i] = reps.get(i);

        return new EquivalenceClasses(classId, repsArr, sizesArr);
    }

    public static final class VertexSignature implements Comparable<VertexSignature> {
        public final int root;                 // helpful for debugging (NOT used for compare)
        public final long[][][] layers;        // layers[d][i][k] = kth entry of ith vector in BFS layer d

        private VertexSignature(int root, long[][][] layers) {
            this.root = root;
            this.layers = layers;
        }

        /** For convenience: how many BFS layers were included (including layer0). */
        public int usedLayers() {
            return layers.length;
        }

        @Override
        public int compareTo(VertexSignature o) {
            if (o == null) return 1;
            int L1 = this.layers.length;
            int L2 = o.layers.length;
            if (L1 != L2) return Integer.compare(L1, L2);

            for (int d = 0; d < L1; d++) {
                long[][] a = this.layers[d];
                long[][] b = o.layers[d];
                if (a.length != b.length) return Integer.compare(a.length, b.length);

                for (int i = 0; i < a.length; i++) {
                    int c = lexCompareLongVec(a[i], b[i]);
                    if (c != 0) return c;
                }
            }
            return 0;
        }

        @Override
        public boolean equals(Object obj) {
            if (!(obj instanceof VertexSignature)) return false;
            return compareTo((VertexSignature) obj) == 0;
        }

        @Override
        public int hashCode() {
            // You said "no hashing" as a matching mechanism. But Java collections need hashCode.
            // We keep it deterministic and consistent with equals, but do not use it for signature identity.
            // (You can still use compareTo/equals for comparisons.)
            int h = 1;
            h = 31 * h + layers.length;
            for (long[][] layer : layers) {
                h = 31 * h + layer.length;
                for (long[] v : layer) {
                    h = 31 * h + Arrays.hashCode(v);
                }
            }
            return h;
        }

        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder("VertexSignature{root=").append(root)
                    .append(", layers=").append(layers.length).append("}");
            return sb.toString();
        }
    }

    public static final class GraphSignature implements Comparable<GraphSignature> {
        public final VertexSignature[] sortedVertexSignatures;

        private GraphSignature(VertexSignature[] sortedVertexSignatures) {
            this.sortedVertexSignatures = Objects.requireNonNull(sortedVertexSignatures, "sortedVertexSignatures null");
        }

        @Override
        public int compareTo(GraphSignature o) {
            if (o == null) return 1;
            int n1 = this.sortedVertexSignatures.length;
            int n2 = o.sortedVertexSignatures.length;
            if (n1 != n2) return Integer.compare(n1, n2);

            for (int i = 0; i < n1; i++) {
                int c = this.sortedVertexSignatures[i].compareTo(o.sortedVertexSignatures[i]);
                if (c != 0) return c;
            }
            return 0;
        }

        @Override
        public boolean equals(Object obj) {
            if (!(obj instanceof GraphSignature)) return false;
            return compareTo((GraphSignature) obj) == 0;
        }

        @Override
        public int hashCode() {
            // Same note as VertexSignature.hashCode()
            return Arrays.hashCode(sortedVertexSignatures);
        }

        @Override
        public String toString() {
            return "GraphSignature{n=" + sortedVertexSignatures.length + "}";
        }
    }

    public static final class EquivalenceClasses {
        public final int[] classIdPerVertex;
        public final int[] representative;
        public final int[] sizes;
        public final int nClasses;

        private EquivalenceClasses(int[] classIdPerVertex, int[] representative, int[] sizes) {
            this.classIdPerVertex = classIdPerVertex;
            this.representative = representative;
            this.sizes = sizes;
            this.nClasses = sizes.length;
        }

        @Override
        public String toString() {

            StringBuilder sb = new StringBuilder();
            sb.append("EquivalenceClasses{nClasses=").append(nClasses).append("}\n");

            // build buckets
            List<List<Integer>> buckets = new ArrayList<>();
            for (int c = 0; c < nClasses; c++) {
                buckets.add(new ArrayList<>(sizes[c]));
            }

            for (int v = 0; v < classIdPerVertex.length; v++) {
                buckets.get(classIdPerVertex[v]).add(v);
            }

            // optional: sort classes by descending size (VERY useful)
            List<Integer> order = new ArrayList<>();
            for (int c = 0; c < nClasses; c++) order.add(c);

            order.sort((a,b) -> Integer.compare(sizes[b], sizes[a]));

            for (int c : order) {
                Collections.sort(buckets.get(c));
                sb.append("  class ")
                  .append(c)
                  .append(" (rep=").append(representative[c])
                  .append(", size=").append(sizes[c])
                  .append("): ")
                  .append(buckets.get(c))
                  .append("\n");
            }

            return sb.toString();
        }
    }

    private static VertexSignature buildVertexSignature(
            Graph<Integer, DefaultEdge> g,
            long[][] qCoeffs,
            int root,
            int bfsLayers
    ) {
        int n = g.vertexSet().size();
        if (root < 0 || root >= n) throw new IllegalArgumentException("root out of range: " + root);

        // BFS distances + per-layer vertex lists
        int maxDepthWanted = bfsLayers - 1; // layer0 is root

        int[] dist = new int[n];
        Arrays.fill(dist, -1);
        dist[root] = 0;

        ArrayDeque<Integer> dq = new ArrayDeque<>();
        dq.add(root);

        // We don’t know diameter; we stop when queue exhausted or reached maxDepthWanted.
        while (!dq.isEmpty()) {
            int u = dq.removeFirst();
            int du = dist[u];
            if (du >= maxDepthWanted) continue;

            for (DefaultEdge e : g.edgesOf(u)) {
                int v = Graphs.getOppositeVertex(g, e, u);
                if (dist[v] == -1) {
                    dist[v] = du + 1;
                    dq.addLast(v);
                }
            }
        }

        int maxDistFound = 0;
        for (int d : dist) if (d > maxDistFound) maxDistFound = d;

        int usedLayers = Math.min(bfsLayers, maxDistFound + 1);

        // Collect vertices by distance 0..usedLayers-1
        @SuppressWarnings("unchecked")
        ArrayList<Integer>[] byLayer = new ArrayList[usedLayers];
        for (int i = 0; i < usedLayers; i++) byLayer[i] = new ArrayList<>();

        for (int v = 0; v < n; v++) {
            int d = dist[v];
            if (d >= 0 && d < usedLayers) byLayer[d].add(v);
        }

        // Convert each layer into sorted list of coeff vectors (layer 0 will have just root).
        long[][][] layers = new long[usedLayers][][];

        for (int d = 0; d < usedLayers; d++) {
            ArrayList<Integer> verts = byLayer[d];

            // Sort vertices in a layer by their coeff vectors lex, tie by vertex id for determinism.
            verts.sort((a, b) -> {
                int c = lexCompareLongVec(qCoeffs[a], qCoeffs[b]);
                if (c != 0) return c;
                return Integer.compare(a, b);
            });

            long[][] vecs = new long[verts.size()][];
            for (int i = 0; i < verts.size(); i++) {
                vecs[i] = qCoeffs[verts.get(i)];
            }
            layers[d] = vecs;
        }

        return new VertexSignature(root, layers);
    }

    /** Lex compare long vectors (by length then entries). */
    public static int lexCompareLongVec(long[] a, long[] b) {
        if (a == b) return 0;
        if (a == null) return (b == null) ? 0 : -1;
        if (b == null) return 1;

        if (a.length != b.length) return Integer.compare(a.length, b.length);
        for (int i = 0; i < a.length; i++) {
            long x = a[i], y = b[i];
            if (x != y) return (x < y) ? -1 : 1;
        }
        return 0;
    }
}
