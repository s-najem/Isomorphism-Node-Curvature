package iso;

import org.jgrapht.Graph;
import org.jgrapht.GraphType;
import org.jgrapht.Graphs;
import org.jgrapht.alg.isomorphism.VF2GraphIsomorphismInspector;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;
import org.jgrapht.graph.builder.GraphTypeBuilder;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

import java.io.BufferedWriter;
import java.nio.charset.StandardCharsets;

public final class GraphUtils 
{
    private GraphUtils() {}

    /* ============================================================
     *  Graph writing (NEW)
     * ============================================================ */

    public enum GraphWriteFormat {
        DIMACS,        // p edge n m; e u v (typically 1-indexed)
        TXT_EDGELIST   // "n m" then "u v" (0-indexed)
    }

    public static void writeGraphToFile(
            Graph<Integer, DefaultEdge> g,
            Path outPath,
            GraphWriteFormat format,
            boolean dimacsOneIndexed
    ) throws IOException {
        requireCanonicalZeroToNMinus1(g);
        Objects.requireNonNull(outPath, "outPath is null");
        Objects.requireNonNull(format, "format is null");

        int n = g.vertexSet().size();

        // Collect edges canonically as (a<b), sort for deterministic output.
        ArrayList<int[]> edges = new ArrayList<>(g.edgeSet().size());
        for (DefaultEdge e : g.edgeSet()) {
            int u = g.getEdgeSource(e);
            int v = g.getEdgeTarget(e);
            if (u == v) continue; // keep simple
            int a = Math.min(u, v);
            int b = Math.max(u, v);
            edges.add(new int[]{a, b});
        }
        edges.sort(java.util.Comparator.<int[]>comparingInt(x -> x[0]).thenComparingInt(x -> x[1]));

        int m = edges.size();

        try (BufferedWriter w = Files.newBufferedWriter(outPath, StandardCharsets.UTF_8)) {
            switch (format) {
                case DIMACS -> {
                    // Optional comment line (safe to remove)
                    w.write("c written by GraphUtils");
                    w.newLine();

                    w.write("p edge " + n + " " + m);
                    w.newLine();

                    int shift = dimacsOneIndexed ? 1 : 0;
                    for (int[] e : edges) {
                        int u = e[0] + shift;
                        int v = e[1] + shift;
                        w.write("e " + u + " " + v);
                        w.newLine();
                    }
                }
                case TXT_EDGELIST -> {
                    // Header: n m
                    w.write(n + " " + m);
                    w.newLine();

                    // 0-indexed edges
                    for (int[] e : edges) {
                        w.write(e[0] + " " + e[1]);
                        w.newLine();
                    }
                }
                default -> throw new IllegalArgumentException("Unknown format: " + format);
            }
        }
    }

    public static void writeDimacs(Path outPath, Graph<Integer, DefaultEdge> g) throws IOException {
        writeGraphToFile(g, outPath, GraphWriteFormat.DIMACS, true);
    }

    public static void writeTxtEdgelist(Path outPath, Graph<Integer, DefaultEdge> g) throws IOException {
        writeGraphToFile(g, outPath, GraphWriteFormat.TXT_EDGELIST, false);
    }

    public static void writeGexf(java.nio.file.Path outPath,
                                 org.jgrapht.Graph<Integer, org.jgrapht.graph.DefaultEdge> g)
            throws java.io.IOException {

        requireCanonicalZeroToNMinus1(g);
        java.util.Objects.requireNonNull(outPath, "outPath is null");

        int n = g.vertexSet().size();

        // Collect edges canonically (a<b) and sort
        java.util.ArrayList<int[]> edges = new java.util.ArrayList<>(g.edgeSet().size());
        for (org.jgrapht.graph.DefaultEdge e : g.edgeSet()) {
            int u = g.getEdgeSource(e);
            int v = g.getEdgeTarget(e);
            if (u == v) continue;

            int a = Math.min(u, v);
            int b = Math.max(u, v);
            edges.add(new int[]{a, b});
        }

        edges.sort(java.util.Comparator
                .<int[]>comparingInt(x -> x[0])
                .thenComparingInt(x -> x[1]));

        try (java.io.BufferedWriter w =
                     java.nio.file.Files.newBufferedWriter(outPath, java.nio.charset.StandardCharsets.UTF_8)) {

            w.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
            w.write("<gexf xmlns=\"http://www.gexf.net/1.2draft\" version=\"1.2\">\n");
            w.write("  <graph mode=\"static\" defaultedgetype=\"undirected\">\n");

            // Nodes
            w.write("    <nodes>\n");
            for (int v = 0; v < n; v++) {
                w.write("      <node id=\"" + v + "\" label=\"" + v + "\" />\n");
            }
            w.write("    </nodes>\n");

            // Edges
            w.write("    <edges>\n");
            int id = 0;
            for (int[] e : edges) {
                w.write("      <edge id=\"" + (id++) +
                        "\" source=\"" + e[0] +
                        "\" target=\"" + e[1] + "\" />\n");
            }
            w.write("    </edges>\n");

            w.write("  </graph>\n");
            w.write("</gexf>\n");
        }
    }

    public enum TwinKind {
        NONE,
        TRUE_TWINS,   // clique internally + same neighbors outside class
        FALSE_TWINS   // independent internally + same neighbors outside class
    }

    public static boolean areTrueTwins(Graph<Integer, DefaultEdge> g, int u, int v, int universeEndExclusive) {
        requireCanonicalZeroToNMinus1(g);
        if (u == v) return false;
        if (universeEndExclusive < 0 || universeEndExclusive > g.vertexSet().size()) {
            throw new IllegalArgumentException("Invalid universeEndExclusive=" + universeEndExclusive);
        }
        if (u < 0 || v < 0 || u >= universeEndExclusive || v >= universeEndExclusive) {
            throw new IllegalArgumentException("u/v out of universe range");
        }

        if (!g.containsEdge(u, v)) return false;

        for (int w = 0; w < universeEndExclusive; w++) {
            if (w == u || w == v) continue;
            boolean auw = g.containsEdge(u, w);
            boolean avw = g.containsEdge(v, w);
            if (auw != avw) return false;
        }
        return true;
    }

    public static boolean areFalseTwins(Graph<Integer, DefaultEdge> g, int u, int v, int universeEndExclusive) {
        requireCanonicalZeroToNMinus1(g);
        if (u == v) return false;
        if (universeEndExclusive < 0 || universeEndExclusive > g.vertexSet().size()) {
            throw new IllegalArgumentException("Invalid universeEndExclusive=" + universeEndExclusive);
        }
        if (u < 0 || v < 0 || u >= universeEndExclusive || v >= universeEndExclusive) {
            throw new IllegalArgumentException("u/v out of universe range");
        }

        if (g.containsEdge(u, v)) return false;

        for (int w = 0; w < universeEndExclusive; w++) {
            if (w == u || w == v) continue;
            boolean auw = g.containsEdge(u, w);
            boolean avw = g.containsEdge(v, w);
            if (auw != avw) return false;
        }
        return true;
    }

    public static TwinKind twinKindOfClass(
            Graph<Integer, DefaultEdge> g,
            List<Integer> vertices,
            int universeEndExclusive
    ) {
        requireCanonicalZeroToNMinus1(g);
        Objects.requireNonNull(vertices, "vertices is null");

        if (universeEndExclusive < 0 || universeEndExclusive > g.vertexSet().size()) {
            throw new IllegalArgumentException("Invalid universeEndExclusive=" + universeEndExclusive);
        }
        if (vertices.size() <= 1) return TwinKind.NONE;

        // Validate vertices and build membership mask (only within universe)
        boolean[] inClass = new boolean[universeEndExclusive];
        for (int x : vertices) {
            if (x < 0 || x >= universeEndExclusive) {
                throw new IllegalArgumentException("Class vertex out of universe: " + x);
            }
            inClass[x] = true;
        }

        // Pick representative
        int rep = vertices.get(0);

        // Record rep adjacency to outside vertices
        boolean[] repAdjOutside = new boolean[universeEndExclusive];
        for (int w = 0; w < universeEndExclusive; w++) {
            if (inClass[w]) continue;
            repAdjOutside[w] = g.containsEdge(rep, w);
        }

        // Check all others match rep adjacency to outside
        for (int i = 1; i < vertices.size(); i++) {
            int x = vertices.get(i);
            for (int w = 0; w < universeEndExclusive; w++) {
                if (inClass[w]) continue;
                boolean axw = g.containsEdge(x, w);
                if (axw != repAdjOutside[w]) return TwinKind.NONE;
            }
        }

        // Check internal structure: clique vs independent
        boolean allEdges = true;
        boolean noEdges = true;

        for (int i = 0; i < vertices.size(); i++) {
            int a = vertices.get(i);
            for (int j = i + 1; j < vertices.size(); j++) {
                int b = vertices.get(j);
                boolean e = g.containsEdge(a, b);
                allEdges &= e;
                noEdges &= !e;
                if (!allEdges && !noEdges) return TwinKind.NONE;
            }
        }

        if (allEdges) return TwinKind.TRUE_TWINS;
        if (noEdges) return TwinKind.FALSE_TWINS;
        return TwinKind.NONE;
    }

    public static void requireCanonicalZeroToNMinus1(Graph<Integer, DefaultEdge> g) {
        Objects.requireNonNull(g, "graph is null");

        int n = g.vertexSet().size();

        // Fast range check + detect duplicates/extras
        for (Integer v : g.vertexSet()) {
            if (v == null) throw new IllegalArgumentException("Graph contains null vertex.");
            if (v < 0 || v >= n) {
                throw new IllegalArgumentException(
                        "Vertex labels must be exactly 0..n-1; found vertex " + v + " with n=" + n
                );
            }
        }

        // Ensure no gaps: each 0..n-1 must exist
        for (int i = 0; i < n; i++) {
            if (!g.containsVertex(i)) {
                throw new IllegalArgumentException(
                        "Vertex labels must be exactly 0..n-1; missing vertex " + i + " with n=" + n
                );
            }
        }
    }

    private static void requireVertexInGraph(Graph<Integer, DefaultEdge> g, int v, String name) {
        if (!g.containsVertex(v)) {
            throw new IllegalArgumentException(name + " not in graph: " + v);
        }
    }

    private static void requirePermutation(int[] perm) {
        Objects.requireNonNull(perm, "perm is null");

        int n = perm.length;
        boolean[] seen = new boolean[n];
        for (int x : perm) {
            if (x < 0 || x >= n) {
                throw new IllegalArgumentException("perm is not a permutation of 0..n-1; out-of-range value: " + x);
            }
            if (seen[x]) {
                throw new IllegalArgumentException("perm is not a permutation of 0..n-1; duplicate value: " + x);
            }
            seen[x] = true;
        }
    }

    public static Graph<Integer, DefaultEdge> readGraphAuto(Path path) throws IOException {
        List<String> lines = readAllNonEmptyLines(path);

        // Try DIMACS/bliss: if any line starts with "p " or "e "
        boolean hasP = false, hasE = false;
        for (String ln : lines) {
            String t = ln.trim();
            if (t.isEmpty()) continue;
            if (t.startsWith("p ")) hasP = true;
            if (t.startsWith("e ")) hasE = true;
        }
        if (hasP || hasE) {
            Graph<Integer, DefaultEdge> g = readDimacsStyle(lines);
            requireCanonicalZeroToNMinus1(g);
            return g;
        }

        // Try header "n m" on first non-comment line, then exactly m edges
        Graph<Integer, DefaultEdge> maybeHeader = tryReadHeaderEdgelist(lines);
        if (maybeHeader != null) {
            requireCanonicalZeroToNMinus1(maybeHeader);
            return maybeHeader;
        }

        // Fallback: flexible edgelist
        Graph<Integer, DefaultEdge> g = readFlexibleEdgelist(lines);
        requireCanonicalZeroToNMinus1(g);
        return g;
    }

    private static List<String> readAllNonEmptyLines(Path path) throws IOException {
        List<String> out = new ArrayList<>();
        try (BufferedReader br = Files.newBufferedReader(path)) {
            String ln;
            while ((ln = br.readLine()) != null) {
                ln = ln.trim();
                if (ln.isEmpty()) continue;
                out.add(ln);
            }
        }
        return out;
    }

    private static boolean isComment(String ln) {
        String t = ln.trim();
        return t.startsWith("#") || t.startsWith("%") || t.startsWith("c") || t.startsWith("C");
    }

    private static Graph<Integer, DefaultEdge> readDimacsStyle(List<String> lines) {
        Integer nDeclared = null;
        List<int[]> edges = new ArrayList<>();

        for (String ln : lines) {
            String t = ln.trim();
            if (t.isEmpty() || isComment(t)) continue;
            if (t.startsWith("p ")) {
                // examples: "p edge n m" or "p sp n m" or "p <something> n m"
                String[] parts = t.split("\\s+");
                // find last two ints
                List<Integer> ints = new ArrayList<>();
                for (String p : parts) {
                    try { ints.add(Integer.parseInt(p)); } catch (NumberFormatException ignored) {}
                }
                if (ints.size() >= 2) {
                    nDeclared = ints.get(ints.size() - 2);
                }
            } else if (t.startsWith("e ")) {
                String[] parts = t.split("\\s+");
                if (parts.length < 3) continue;
                int u = Integer.parseInt(parts[1]);
                int v = Integer.parseInt(parts[2]);
                edges.add(new int[]{u, v});
            } else {
                // ignore other record types
            }
        }

        // DIMACS is usually 1-indexed. Apply shift heuristic if min==1 and no 0 present.
        return buildGraphFromEdges(edges, nDeclared, true);
    }

    private static Graph<Integer, DefaultEdge> tryReadHeaderEdgelist(List<String> lines) {
        int idx = 0;
        while (idx < lines.size() && isComment(lines.get(idx))) idx++;
        if (idx >= lines.size()) return null;

        String first = lines.get(idx).trim();
        String[] parts = first.split("\\s+");
        if (parts.length != 2) return null;

        int n, m;
        try {
            n = Integer.parseInt(parts[0]);
            m = Integer.parseInt(parts[1]);
        } catch (NumberFormatException e) {
            return null;
        }
        if (n <= 0 || m < 0) return null;

        List<int[]> edges = new ArrayList<>();
        for (int j = idx + 1; j < lines.size() && edges.size() < m; j++) {
            String ln = lines.get(j).trim();
            if (ln.isEmpty() || isComment(ln)) continue;
            String[] p = ln.split("\\s+");
            if (p.length < 2) continue;
            try {
                int u = Integer.parseInt(p[0]);
                int v = Integer.parseInt(p[1]);
                edges.add(new int[]{u, v});
            } catch (NumberFormatException ignored) {}
        }
        if (edges.size() < m) return null;

        // Important: if this "header format" guess is wrong, DO NOT hard-fail; fall back.
        try {
            return buildGraphFromEdges(edges, n, true);
        } catch (IllegalArgumentException ex) {
            return null;
        }
    }

    private static Graph<Integer, DefaultEdge> readFlexibleEdgelist(List<String> lines) {
        List<int[]> edges = new ArrayList<>();
        for (String ln : lines) {
            String t = ln.trim();
            if (t.isEmpty() || isComment(t)) continue;
            String[] p = t.split("\\s+");
            if (p.length < 2) continue;
            try {
                int u = Integer.parseInt(p[0]);
                int v = Integer.parseInt(p[1]);
                edges.add(new int[]{u, v});
            } catch (NumberFormatException ignored) {}
        }
        // n unknown; infer.
        return buildGraphFromEdges(edges, null, true);
    }

    private static Graph<Integer, DefaultEdge> buildGraphFromEdges(List<int[]> edges, Integer nDeclared, boolean applyShiftHeuristic) {
        if (edges == null) throw new IllegalArgumentException("edges is null");

        if (edges.isEmpty()) {
            if (nDeclared == null || nDeclared <= 0) {
                throw new IllegalArgumentException("No edges parsed and n is unknown/invalid.");
            }
            Graph<Integer, DefaultEdge> g = new SimpleGraph<>(DefaultEdge.class);
            for (int i = 0; i < nDeclared; i++) g.addVertex(i);
            return g;
        }

        // collect node labels
        IntSummaryStatistics stats = new IntSummaryStatistics();
        Set<Integer> nodes = new HashSet<>();
        boolean hasZero = false;
        for (int[] e : edges) {
            if (e == null || e.length < 2) continue;
            nodes.add(e[0]); nodes.add(e[1]);
            stats.accept(e[0]); stats.accept(e[1]);
            if (e[0] == 0 || e[1] == 0) hasZero = true;
        }

        int min = stats.getMin();
        int max = stats.getMax();

        boolean shift = applyShiftHeuristic && (min == 1) && !hasZero;
        if (shift) {
            List<int[]> shifted = new ArrayList<>(edges.size());
            for (int[] e : edges) shifted.add(new int[]{e[0] - 1, e[1] - 1});
            edges = shifted;

            // recompute nodes
            nodes.clear();
            max = Integer.MIN_VALUE;
            for (int[] e : edges) {
                nodes.add(e[0]); nodes.add(e[1]);
                max = Math.max(max, Math.max(e[0], e[1]));
            }
            min = 0;
        }

        int nInferred = (nDeclared != null) ? nDeclared : (max + 1);

        if (min != 0) {
            throw new IllegalArgumentException("Vertex labels must start at 0 after shift heuristic. min=" + min);
        }
        if (nInferred <= 0) throw new IllegalArgumentException("Invalid inferred n=" + nInferred);

        // If nDeclared provided, require that all labels are in [0, n-1]
        for (int v : nodes) {
            if (v < 0 || v >= nInferred) {
                throw new IllegalArgumentException("Vertex label out of range: " + v + " for n=" + nInferred);
            }
        }

        Graph<Integer, DefaultEdge> g = new SimpleGraph<>(DefaultEdge.class);
        for (int i = 0; i < nInferred; i++) g.addVertex(i);

        // Add edges; enforce simple (no self loops)
        for (int[] e : edges) {
            int u = e[0], v = e[1];
            if (u == v) continue; // stay simple
            g.addEdge(u, v);      // SimpleGraph prevents multiedges
        }

        return g;
    }

    public static Graph<Integer, DefaultEdge> permute(Graph<Integer, DefaultEdge> g, int[] perm) {
        requireCanonicalZeroToNMinus1(g);
        if (perm.length != g.vertexSet().size()) {
            throw new IllegalArgumentException("perm length mismatch: perm.length=" + perm.length +
                    " n=" + g.vertexSet().size());
        }
        requirePermutation(perm);

        int n = g.vertexSet().size();

        Graph<Integer, DefaultEdge> gp = new SimpleGraph<>(DefaultEdge.class);
        for (int i = 0; i < n; i++) gp.addVertex(i);

        for (DefaultEdge e : g.edgeSet()) {
            int u = g.getEdgeSource(e);
            int v = g.getEdgeTarget(e);
            int up = perm[u];
            int vp = perm[v];
            if (up == vp) continue;
            gp.addEdge(up, vp);
        }
        return gp;
    }

    public static int[] randomPermutation(int n, long seed) {
        if (n < 0) throw new IllegalArgumentException("n must be >= 0");
        int[] p = new int[n];
        for (int i = 0; i < n; i++) p[i] = i;
        Random rnd = new Random(seed);
        for (int i = n - 1; i > 0; i--) {
            int j = rnd.nextInt(i + 1);
            int tmp = p[i]; p[i] = p[j]; p[j] = tmp;
        }
        return p;
    }

    public static Graph<Integer, DefaultEdge> attachClique(
            Graph<Integer, DefaultEdge> g,
            int attachVertex,
            int cliqueSize,
            boolean connectAttachToAll
    ) {
        requireCanonicalZeroToNMinus1(g);
        if (cliqueSize < 0) throw new IllegalArgumentException("cliqueSize must be >= 0");
        requireVertexInGraph(g, attachVertex, "attachVertex");

        int n = g.vertexSet().size();

        Graph<Integer, DefaultEdge> out = copyCanonicalGraph(g);

        // Add new clique vertices
        int[] newVerts = new int[cliqueSize];
        for (int i = 0; i < cliqueSize; i++) {
            int vNew = n + i;
            newVerts[i] = vNew;
            out.addVertex(vNew);
        }

        // Add clique internal edges
        for (int i = 0; i < cliqueSize; i++) {
            for (int j = i + 1; j < cliqueSize; j++) {
                out.addEdge(newVerts[i], newVerts[j]);
            }
        }

        // Connect clique to attachVertex
        if (connectAttachToAll) {
            for (int vNew : newVerts) out.addEdge(attachVertex, vNew);
        } else if (cliqueSize > 0) {
            out.addEdge(attachVertex, newVerts[0]);
        }

        return out;
    }

    public static Graph<Integer, DefaultEdge> attachCliqueAsSuperClique(
            Graph<Integer, DefaultEdge> g,
            int attachVertex,
            int cliqueSize
    ) {
        return attachClique(g, attachVertex, cliqueSize, true);
    }

    public static Graph<Integer, DefaultEdge> addUniversalVertex(Graph<Integer, DefaultEdge> g) {
        requireCanonicalZeroToNMinus1(g);

        int n = g.vertexSet().size();

        Graph<Integer, DefaultEdge> out = copyCanonicalGraph(g);

        int uNew = n;
        out.addVertex(uNew);

        for (int v = 0; v < n; v++) out.addEdge(uNew, v);

        return out;
    }

    public static Graph<Integer, DefaultEdge> addUniversalPath(Graph<Integer, DefaultEdge> g, int k) {
        requireCanonicalZeroToNMinus1(g);
        if (k < 0) throw new IllegalArgumentException("k must be >= 0");

        int n = g.vertexSet().size();
        Graph<Integer, DefaultEdge> out = copyCanonicalGraph(g);

        // Add new vertices
        for (int i = 0; i < k; i++) {
            out.addVertex(n + i);
        }

        // Make each new vertex universal to the ORIGINAL vertices (0..n-1)
        for (int i = 0; i < k; i++) {
            int u = n + i;
            for (int v = 0; v < n; v++) {
                out.addEdge(u, v);
            }
        }

        // Connect new vertices as a path: (n)-(n+1)-...-(n+k-1)
        for (int i = 0; i + 1 < k; i++) {
            out.addEdge(n + i, n + i + 1);
        }

        return out;
    }

    public static Graph<Integer, DefaultEdge> attachPath(
            Graph<Integer, DefaultEdge> g,
            int attachVertex,
            int pathSize
    ) {
        requireCanonicalZeroToNMinus1(g);
        if (pathSize < 0) throw new IllegalArgumentException("pathSize must be >= 0");
        requireVertexInGraph(g, attachVertex, "attachVertex");

        int n = g.vertexSet().size();

        Graph<Integer, DefaultEdge> out = copyCanonicalGraph(g);

        int[] p = new int[pathSize];
        for (int i = 0; i < pathSize; i++) {
            int vNew = n + i;
            p[i] = vNew;
            out.addVertex(vNew);
        }

        if (pathSize > 0) {
            out.addEdge(attachVertex, p[0]);
            for (int i = 0; i + 1 < pathSize; i++) {
                out.addEdge(p[i], p[i + 1]);
            }
        }

        return out;
    }

    public static Graph<Integer, DefaultEdge> attachCliqueAndPathAsSuperClique(
            Graph<Integer, DefaultEdge> g,
            int attachVertex,
            int K
    ) {
        Graph<Integer, DefaultEdge> withClique = attachCliqueAsSuperClique(g, attachVertex, K);
        return attachPath(withClique, attachVertex, K);
    }

    public static Graph<Integer, DefaultEdge> attachCliqueAndTwoPathsAsSuperClique(
            Graph<Integer, DefaultEdge> g,
            int attachVertex,
            int K
    ) {
        Graph<Integer, DefaultEdge> withClique = attachCliqueAsSuperClique(g, attachVertex, K);
        Graph<Integer, DefaultEdge> withPath = attachPath(withClique, attachVertex, K);
        return attachPath(withPath, attachVertex, K + 1);
    }

    private static Graph<Integer, DefaultEdge> copyCanonicalGraph(Graph<Integer, DefaultEdge> g) {
        // g already validated by requireCanonicalZeroToNMinus1 in callers
        int n = g.vertexSet().size();

        Graph<Integer, DefaultEdge> out = new SimpleGraph<>(DefaultEdge.class);
        for (int v = 0; v < n; v++) out.addVertex(v);

        for (DefaultEdge e : g.edgeSet()) {
            int u = g.getEdgeSource(e);
            int v = g.getEdgeTarget(e);
            if (u == v) continue; // keep simple
            out.addEdge(u, v);
        }
        return out;
    }

    public static Graph<Integer, DefaultEdge> subdivideEveryEdge(Graph<Integer, DefaultEdge> g, int k) {
        requireCanonicalZeroToNMinus1(g);
        if (k < 0) throw new IllegalArgumentException("k must be >= 0");

        int n = g.vertexSet().size();

        // Start output with original vertices
        Graph<Integer, DefaultEdge> out = new SimpleGraph<>(DefaultEdge.class);
        for (int v = 0; v < n; v++) out.addVertex(v);

        // If k==0, just copy edges (but still deterministically)
        if (k == 0) {
            // Collect edges canonically (a<b)
            ArrayList<int[]> edges = new ArrayList<>(g.edgeSet().size());
            for (DefaultEdge e : g.edgeSet()) {
                int u = g.getEdgeSource(e);
                int v = g.getEdgeTarget(e);
                if (u == v) continue;
                int a = Math.min(u, v);
                int b = Math.max(u, v);
                edges.add(new int[]{a, b});
            }
            edges.sort(java.util.Comparator.<int[]>comparingInt(x -> x[0]).thenComparingInt(x -> x[1]));

            for (int[] edge : edges) {
                out.addEdge(edge[0], edge[1]);
            }
            return out;
        }

        int nextNew = n;

        // Collect original edges canonically (a<b) and sort to ensure deterministic subdivision order.
        ArrayList<int[]> edges = new ArrayList<>(g.edgeSet().size());
        for (DefaultEdge e : g.edgeSet()) {
            int u = g.getEdgeSource(e);
            int v = g.getEdgeTarget(e);

            // SimpleGraph should not have self-loops, but keep it safe.
            if (u == v) continue;

            int a = Math.min(u, v);
            int b = Math.max(u, v);
            edges.add(new int[]{a, b});
        }

        edges.sort(java.util.Comparator.<int[]>comparingInt(x -> x[0]).thenComparingInt(x -> x[1]));

        // Subdivide edges in deterministic order.
        for (int[] edge : edges) {
            int u = edge[0];
            int v = edge[1];

            int prev = u;

            // Insert k new vertices on the edge u--v: u--x1--x2--...--xk--v
            for (int i = 0; i < k; i++) {
                int x = nextNew++;
                out.addVertex(x);
                out.addEdge(prev, x);
                prev = x;
            }

            out.addEdge(prev, v);
        }

        return out;
    }
    
    public static Graph<Integer, DefaultEdge> relabelBFS(Graph<Integer, DefaultEdge> g, int root) {
        requireCanonicalZeroToNMinus1(g);

        int n = g.vertexSet().size();

        Map<Integer, Integer> newId = new HashMap<>(n);
        Queue<Integer> q = new ArrayDeque<>();

        q.add(root);
        newId.put(root, 0);

        int next = 1;

        while (!q.isEmpty()) {
            int u = q.poll();

            List<Integer> nbrs = Graphs.neighborListOf(g, u);
            Collections.sort(nbrs); // deterministic traversal

            for (int v : nbrs) {
                if (!newId.containsKey(v)) {
                    newId.put(v, next++);
                    q.add(v);
                }
            }
        }

        // rebuild graph
        Graph<Integer, DefaultEdge> out = new SimpleGraph<>(DefaultEdge.class);
        for (int i = 0; i < n; i++) out.addVertex(i);

        for (DefaultEdge e : g.edgeSet()) {
            int u = g.getEdgeSource(e);
            int v = g.getEdgeTarget(e);
            out.addEdge(newId.get(u), newId.get(v));
        }

        return out;
    }
    
    public static boolean areIsomorphic(Graph<Integer, DefaultEdge> g1,
            Graph<Integer, DefaultEdge> g2) {
    	// If you have vertex/edge labels you care about, see section (2) below.
    	VF2GraphIsomorphismInspector<Integer, DefaultEdge> insp = new VF2GraphIsomorphismInspector<>(g1, g2);
    	return insp.isomorphismExists();
    }

    public static Graph<Integer, DefaultEdge> attachCycleWithPendants(
            Graph<Integer, DefaultEdge> g,
            int attachVertex,
            int cycleSize
    ) {
        requireCanonicalZeroToNMinus1(g);
        requireVertexInGraph(g, attachVertex, "attachVertex");
        if (cycleSize < 0) throw new IllegalArgumentException("cycleSize must be >= 0");
        if (cycleSize > 0 && cycleSize < 3) {
            throw new IllegalArgumentException("cycleSize must be 0 or >= 3 (a cycle needs at least 3 vertices)");
        }

        int n = g.vertexSet().size();
        Graph<Integer, DefaultEdge> out = copyCanonicalGraph(g);

        if (cycleSize == 0) return out; // attach nothing

        int[] cycle = new int[cycleSize];
        int[] pendant = new int[cycleSize];

        // Add cycle vertices
        for (int i = 0; i < cycleSize; i++) {
            cycle[i] = n + i;
            out.addVertex(cycle[i]);
        }

        // Add pendant vertices
        for (int i = 0; i < cycleSize; i++) {
            pendant[i] = n + cycleSize + i;
            out.addVertex(pendant[i]);
        }

        // Connect attachment vertex to the first cycle vertex
        out.addEdge(attachVertex, cycle[0]);

        // Add cycle edges
        for (int i = 0; i < cycleSize; i++) {
            int a = cycle[i];
            int b = cycle[(i + 1) % cycleSize];
            out.addEdge(a, b);
        }

        // Add one pendant edge per cycle vertex
        for (int i = 0; i < cycleSize; i++) {
            out.addEdge(cycle[i], pendant[i]);
        }

        return out;
    }
    
    public static Graph<Integer, DefaultEdge> attachCycle(
            Graph<Integer, DefaultEdge> g,
            int attachVertex,
            int cycleSize
    ) {
        requireCanonicalZeroToNMinus1(g);
        requireVertexInGraph(g, attachVertex, "attachVertex");
        if (cycleSize < 0) throw new IllegalArgumentException("cycleSize must be >= 0");
        if (cycleSize > 0 && cycleSize < 3) {
            throw new IllegalArgumentException("cycleSize must be 0 or >= 3 (a cycle needs at least 3 vertices)");
        }

        int n = g.vertexSet().size();
        Graph<Integer, DefaultEdge> out = copyCanonicalGraph(g);

        if (cycleSize == 0) return out; // attach nothing

        int[] cycle = new int[cycleSize];

        // Add cycle vertices
        for (int i = 0; i < cycleSize; i++) {
            cycle[i] = n + i;
            out.addVertex(cycle[i]);
        }

        // Connect attachment vertex to the first cycle vertex
        out.addEdge(attachVertex, cycle[0]);

        // Add cycle edges
        for (int i = 0; i < cycleSize; i++) {
            int a = cycle[i];
            int b = cycle[(i + 1) % cycleSize];
            out.addEdge(a, b);
        }

        return out;
    }

    public static Graph<Integer, DefaultEdge> attachFingerprintLollipop(
            Graph<Integer, DefaultEdge> g,
            int attachVertex,
            int idx,
            int baseTail,
            int baseCycle
    ) {
        requireCanonicalZeroToNMinus1(g);
        requireVertexInGraph(g, attachVertex, "attachVertex");
        if (idx < 0) throw new IllegalArgumentException("idx must be >= 0");
        if (baseTail < 0) throw new IllegalArgumentException("baseTail must be >= 0");
        if (baseCycle < 3) throw new IllegalArgumentException("baseCycle must be >= 3");

        // Make parameters unique per idx
        int tailLen = baseTail + idx;                 // all distinct
        int cycLen  = nextPrime(baseCycle + 2 * idx); // all distinct, primes reduce collisions

        Graph<Integer, DefaultEdge> out = copyCanonicalGraph(g);
        int nextNew = out.vertexSet().size();

        // Tail vertices
        int[] tail = new int[tailLen];
        for (int i = 0; i < tailLen; i++) {
            tail[i] = nextNew++;
            out.addVertex(tail[i]);
        }

        // Cycle vertices
        int[] cyc = new int[cycLen];
        for (int j = 0; j < cycLen; j++) {
            cyc[j] = nextNew++;
            out.addVertex(cyc[j]);
        }

        // Attach: attachVertex -- tailHead (or attachVertex -- cyc[0] if no tail)
        if (tailLen > 0) {
            out.addEdge(attachVertex, tail[0]);
            for (int i = 0; i + 1 < tailLen; i++) out.addEdge(tail[i], tail[i + 1]);
            out.addEdge(tail[tailLen - 1], cyc[0]);
        } else {
            out.addEdge(attachVertex, cyc[0]);
        }

        // Cycle edges
        for (int j = 0; j < cycLen; j++) {
            out.addEdge(cyc[j], cyc[(j + 1) % cycLen]);
        }

        return out;
    }

    private static int nextPrime(int x) {
        if (x <= 2) return 2;
        int p = (x % 2 == 0) ? x + 1 : x;
        while (!isPrime(p)) p += 2;
        return p;
    }

    private static boolean isPrime(int x) {
        if (x < 2) return false;
        if (x % 2 == 0) return x == 2;
        for (int d = 3; (long) d * d <= x; d += 2) {
            if (x % d == 0) return false;
        }
        return true;
    }

    public static Graph<Integer, DefaultEdge> attachSharedClique(
            Graph<Integer, DefaultEdge> g,
            List<Integer> attachVertices,
            int cliqueSize
    ) {
        requireCanonicalZeroToNMinus1(g);
        Objects.requireNonNull(attachVertices, "attachVertices is null");

        if (cliqueSize < 0)
            throw new IllegalArgumentException("cliqueSize must be >= 0");

        if (attachVertices.isEmpty())
            throw new IllegalArgumentException("attachVertices cannot be empty");

        // validate vertices
        boolean[] seen = new boolean[g.vertexSet().size()];
        for (int v : attachVertices) {
            requireVertexInGraph(g, v, "attachVertex");
            if (seen[v])
                throw new IllegalArgumentException("Duplicate vertex in attachVertices: " + v);
            seen[v] = true;
        }

        int n = g.vertexSet().size();
        Graph<Integer, DefaultEdge> out = copyCanonicalGraph(g);

        // create clique vertices
        int[] clique = new int[cliqueSize];
        for (int i = 0; i < cliqueSize; i++) {
            clique[i] = n + i;
            out.addVertex(clique[i]);
        }

        // internal clique edges
        for (int i = 0; i < cliqueSize; i++) {
            for (int j = i + 1; j < cliqueSize; j++) {
                out.addEdge(clique[i], clique[j]);
            }
        }

        // connect clique to all attachment vertices
        for (int vAttach : attachVertices) {
            for (int c : clique) {
                out.addEdge(vAttach, c);
            }
        }

        return out;
    }
    
    public static Graph<Integer, DefaultEdge> attachSharedClique(
            Graph<Integer, DefaultEdge> g,
            int cliqueSize,
            int... attachVertices
    ) {
        List<Integer> list = new ArrayList<>(attachVertices.length);
        for (int v : attachVertices) list.add(v);
        return attachSharedClique(g, list, cliqueSize);
    }
    
    public static Graph<Integer, DefaultEdge> attachSharedCliqueWithDistinctMarkers(
            Graph<Integer, DefaultEdge> g,
            List<Integer> vertices,
            int sharedCliqueSize,
            int basePathLen // e.g. 2
    ) {
        requireCanonicalZeroToNMinus1(g);
        Objects.requireNonNull(vertices, "vertices is null");
        if (vertices.isEmpty()) throw new IllegalArgumentException("vertices cannot be empty");
        if (sharedCliqueSize < 0) throw new IllegalArgumentException("sharedCliqueSize must be >= 0");
        if (basePathLen < 0) throw new IllegalArgumentException("basePathLen must be >= 0");

        // validate + distinct
        boolean[] seen = new boolean[g.vertexSet().size()];
        for (int v : vertices) {
            requireVertexInGraph(g, v, "vertex");
            if (seen[v]) throw new IllegalArgumentException("Duplicate vertex in vertices: " + v);
            seen[v] = true;
        }

        // Step 1: attach ONE clique shared by ALL vertices in the list
        Graph<Integer, DefaultEdge> out = attachSharedClique(g, vertices, sharedCliqueSize);

        // Step 2: attach distinct marker paths to each vertex
        // marker length for vertices[i] = basePathLen + i
        for (int i = 0; i < vertices.size(); i++) {
            int v = vertices.get(i);
            out = attachPath(out, v, basePathLen + i);
        }

        return out;
    }
    
    public static void printNeighbors(String label, Graph<Integer, DefaultEdge> g, int... verts) {
        System.out.println("==== Neighbors in " + label + " ====");
        for (int v : verts) {
            if (!g.containsVertex(v)) {
                System.out.println(label + " v=" + v + " : (vertex not present)");
                continue;
            }
            List<Integer> nbrs = new ArrayList<>();
            for (DefaultEdge e : g.edgesOf(v)) {
                Integer a = g.getEdgeSource(e);
                Integer b = g.getEdgeTarget(e);
                nbrs.add(a.equals(v) ? b : a);
            }
            nbrs.sort(Integer::compareTo);
            System.out.println(label + " v=" + v + " deg=" + nbrs.size() + " nbrs=" + nbrs);
        }
        System.out.println();
    }
    
    public static void printNeighborsOriginalOnly(String label, Graph<Integer, DefaultEdge> g, int nOriginal, int... verts) {
        System.out.println("==== Neighbors in " + label + " (original-only) ====");
        for (int v : verts) {
            if (!g.containsVertex(v)) {
                System.out.println(label + " v=" + v + " : (vertex not present)");
                continue;
            }
            List<Integer> nbrs = new ArrayList<>();
            for (DefaultEdge e : g.edgesOf(v)) {
                Integer a = g.getEdgeSource(e);
                Integer b = g.getEdgeTarget(e);
                int w = a.equals(v) ? b : a;
                if (w >= 0 && w < nOriginal) nbrs.add(w);
            }
            nbrs.sort(Integer::compareTo);
            System.out.println(label + " v=" + v + " deg(original)=" + nbrs.size() + " nbrs=" + nbrs);
        }
        System.out.println();
    }

    public static Graph<Integer, DefaultEdge> copyGraphSameType(Graph<Integer, DefaultEdge> g) {
        GraphType t = g.getType();

        Graph<Integer, DefaultEdge> h =
                GraphTypeBuilder
                        .<Integer, DefaultEdge>forGraphType(t)
                        .edgeClass(DefaultEdge.class)
                        .buildGraph();

        for (Integer v : g.vertexSet()) h.addVertex(v);

        boolean weighted = t.isWeighted();
        for (DefaultEdge e : g.edgeSet()) {
            Integer s = g.getEdgeSource(e);
            Integer tt = g.getEdgeTarget(e);
            DefaultEdge e2 = h.addEdge(s, tt);
            if (weighted && e2 != null) {
                h.setEdgeWeight(e2, g.getEdgeWeight(e));
            }
        }
        return h;
    }
}

