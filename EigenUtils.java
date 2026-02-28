package iso;

import org.ejml.data.Complex_F64;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.factory.DecompositionFactory_DDRM;
import org.ejml.interfaces.decomposition.EigenDecomposition_F64;
import org.jgrapht.Graph;
import org.jgrapht.Graphs;
import org.jgrapht.graph.DefaultEdge;

import java.util.ArrayDeque;
import java.util.Arrays;
import java.util.Deque;
import java.util.Objects;

public final class EigenUtils 
{    
    private EigenUtils() {}

    /** Which symmetric matrix to build from a graph. */
    public enum GraphMatrixType {
        ADJACENCY,           // A
        LAPLACIAN,           // L = D - A
        NORMALIZED_LAPLACIAN // Lsym = I - D^{-1/2} A D^{-1/2}
    }

    public static boolean sameSpectrumQuantized(double[] e1, double[] e2, double scale) {
        if (e1.length != e2.length) return false;

        long[] q1 = EigenUtils.quantize(e1, scale);
        long[] q2 = EigenUtils.quantize(e2, scale);

        Arrays.sort(q1);
        Arrays.sort(q2);

        return Arrays.equals(q1, q2);
    }

    public static boolean nearlyEqual(double a, double b, double eps) {
        double diff = Math.abs(a - b);

        if (a == b) return true; // handles infinities

        double largest = Math.max(Math.abs(a), Math.abs(b));

        return diff <= largest * eps || diff <= eps;
    }

    static int countUniqueSorted(long[] a) {
        if (a.length == 0) return 0;
        int count = 1;
        long prev = a[0];
        for (int i = 1; i < a.length; i++) {
            long x = a[i];
            if (x != prev) {
                count++;
                prev = x;
            }
        }
        return count;
    }

    static int countUnique(long[] a) {
        if (a.length == 0) return 0;
        long[] b = a.clone();
        Arrays.sort(b);
        return countUniqueSorted(b);
    }
    
    public static long[] quantize(double[] evals) {
        return quantize(evals, Constants.EIGEN_DEFAULT_SCALE);
    }

    public static long[] quantize(double[] evals, double scale) {
        if (scale <= 0) {
            throw new IllegalArgumentException("scale must be > 0");
        }

        long[] q = new long[evals.length];
        double clamp = 0.5 / scale;

        for (int i = 0; i < evals.length; i++) {
            double x = evals[i];

            // kill tiny FP noise and negative zero
            if (Math.abs(x) < clamp) {
                q[i] = 0;
                continue;
            }

            double scaled = x * scale;

            // extremely defensive — practically never triggered
            if (scaled > Long.MAX_VALUE || scaled < Long.MIN_VALUE) {
                throw new ArithmeticException("Quantization overflow: " + x);
            }

            q[i] = Math.round(scaled);
        }

        return q;
    }

    public static double[] canonicalize(double[] evals) {
        double[] out = evals.clone();

        for (int i = 0; i < out.length; i++) {
            if (Math.abs(out[i]) < Constants.EIGEN_ZERO_CLAMP) {
                out[i] = 0.0; // removes -0.0
            }
        }

        return out;
    }
    
    public static double[] computeSpectrumFromGraph(
            Graph<Integer, DefaultEdge> g,
            GraphMatrixType type,
            DMatrixRMaj eigenvectorsColsOut
    ) {
        Objects.requireNonNull(g, "graph is null");
        Objects.requireNonNull(type, "type is null");
        Objects.requireNonNull(eigenvectorsColsOut, "eigenvectorsColsOut is null");

        requireSimpleConnectedCanonical(g);
        DMatrixRMaj M = buildMatrixFromGraph(g, type);
        return solveSymmetric(M, eigenvectorsColsOut);
    }

    /**
     * Convenience overload that allocates eigenvector matrix.
     * @return Object[]{ double[] eigenvaluesAscending, DMatrixRMaj eigenvectorsCols }
     */
    public static Object[] computeSpectrumFromGraph(Graph<Integer, DefaultEdge> g, GraphMatrixType type) {
        Objects.requireNonNull(g, "graph is null");
        Objects.requireNonNull(type, "type is null");
        int n = g.vertexSet().size();
        DMatrixRMaj V = new DMatrixRMaj(n, n);
        double[] evals = computeSpectrumFromGraph(g, type, V);
        return new Object[]{ evals, V };
    }

    /** Builds the chosen symmetric matrix (dense) from the graph. Graph must be canonical-labeled. */
    public static DMatrixRMaj buildMatrixFromGraph(Graph<Integer, DefaultEdge> g, GraphMatrixType type) {
        Objects.requireNonNull(g, "graph is null");
        Objects.requireNonNull(type, "type is null");

        requireCanonicalZeroToNMinus1(g); // needed for matrix indexing

        int n = g.vertexSet().size();
        DMatrixRMaj M = new DMatrixRMaj(n, n);

        // degree array (needed for L and Lsym)
        int[] deg = new int[n];
        for (int v = 0; v < n; v++) deg[v] = g.degreeOf(v);

        switch (type) {
            case ADJACENCY -> {
                for (DefaultEdge e : g.edgeSet()) {
                    int u = g.getEdgeSource(e);
                    int v = g.getEdgeTarget(e);
                    if (u == v) continue;
                    M.set(u, v, 1.0);
                    M.set(v, u, 1.0);
                }
                return M;
            }
            case LAPLACIAN -> {
                for (int i = 0; i < n; i++) M.set(i, i, (double) deg[i]);
                for (DefaultEdge e : g.edgeSet()) {
                    int u = g.getEdgeSource(e);
                    int v = g.getEdgeTarget(e);
                    if (u == v) continue;
                    M.set(u, v, -1.0);
                    M.set(v, u, -1.0);
                }
                return M;
            }
            case NORMALIZED_LAPLACIAN -> {
                // For connected simple graphs, degrees should be > 0 everywhere.
                for (int i = 0; i < n; i++) {
                    if (deg[i] <= 0) {
                        throw new IllegalArgumentException(
                                "Normalized Laplacian requires deg(v)>0 for all v; found deg(" + i + ")=" + deg[i]
                        );
                    }
                }

                // Start with identity
                for (int i = 0; i < n; i++) M.set(i, i, 1.0);

                for (DefaultEdge e : g.edgeSet()) {
                    int u = g.getEdgeSource(e);
                    int v = g.getEdgeTarget(e);
                    if (u == v) continue;

                    // Compute once; assign both directions from the same value to guarantee exact symmetry bitwise.
                    double w = 1.0 / Math.sqrt((double) deg[u] * (double) deg[v]);
                    M.set(u, v, -w);
                    M.set(v, u, -w);
                }
                return M;
            }
            default -> throw new IllegalStateException("Unhandled type: " + type);
        }
    }

    /** Run a strong set of checks for your setting: canonical labels, no self-loops, connected. */
    public static void requireSimpleConnectedCanonical(Graph<Integer, DefaultEdge> g) {
        requireCanonicalZeroToNMinus1(g);
        requireNoSelfLoops(g);
        requireConnected(g);
        // Optional quick check: connected simple graph implies m >= n-1
        requireEdgeCountLowerBound(g);
    }

    /** Ensures vertices are exactly {0..n-1} where n = vertexSet().size(). */
    public static void requireCanonicalZeroToNMinus1(Graph<Integer, DefaultEdge> g) {
        Objects.requireNonNull(g, "graph is null");
        int n = g.vertexSet().size();

        for (Integer v : g.vertexSet()) {
            if (v == null) throw new IllegalArgumentException("Graph contains a null vertex.");
            if (v < 0 || v >= n) {
                throw new IllegalArgumentException("Vertex labels must be exactly 0..n-1; found " + v + " with n=" + n);
            }
        }
        for (int i = 0; i < n; i++) {
            if (!g.containsVertex(i)) {
                throw new IllegalArgumentException("Vertex labels must be exactly 0..n-1; missing vertex " + i + " with n=" + n);
            }
        }
    }

    /** Ensures there are no self-loop edges (v,v). */
    public static void requireNoSelfLoops(Graph<Integer, DefaultEdge> g) {
        Objects.requireNonNull(g, "graph is null");
        for (DefaultEdge e : g.edgeSet()) {
            Integer u = g.getEdgeSource(e);
            Integer v = g.getEdgeTarget(e);
            if (u == null || v == null) throw new IllegalArgumentException("Edge has null endpoint: " + e);
            if (u.equals(v)) throw new IllegalArgumentException("Self-loop detected at vertex " + u);
        }
    }

    /** Ensures the graph is connected (treating it as undirected). Requires canonical labels so we can use boolean[]. */
    public static void requireConnected(Graph<Integer, DefaultEdge> g) {
        Objects.requireNonNull(g, "graph is null");
        int n = g.vertexSet().size();
        if (n == 0) throw new IllegalArgumentException("Graph must be non-empty.");
        if (n == 1) return;

        boolean[] seen = new boolean[n];
        Deque<Integer> dq = new ArrayDeque<>();
        dq.add(0);
        seen[0] = true;

        while (!dq.isEmpty()) {
            int u = dq.removeFirst();
            for (DefaultEdge e : g.edgesOf(u)) {
                int v = Graphs.getOppositeVertex(g, e, u);
                if (v < 0 || v >= n) {
                    throw new IllegalArgumentException("Found neighbor with label out of range: " + v);
                }
                if (!seen[v]) {
                    seen[v] = true;
                    dq.addLast(v);
                }
            }
        }

        for (int i = 0; i < n; i++) {
            if (!seen[i]) throw new IllegalArgumentException("Graph is not connected; unreachable vertex: " + i);
        }
    }

    /** Quick necessary condition: connected simple graph on n vertices must satisfy m >= n-1. */
    public static void requireEdgeCountLowerBound(Graph<Integer, DefaultEdge> g) {
        Objects.requireNonNull(g, "graph is null");
        int n = g.vertexSet().size();
        int m = g.edgeSet().size();
        if (n > 0 && m < n - 1) {
            throw new IllegalArgumentException("Too few edges for connected graph: n=" + n + " m=" + m);
        }
    }

    public static void requireConnectedLaplacianSignature(double[] laplacianEvalsAscending, double tol) {
        Objects.requireNonNull(laplacianEvalsAscending, "laplacianEvalsAscending is null");
        if (tol < 0) throw new IllegalArgumentException("tol must be >= 0");
        if (laplacianEvalsAscending.length == 0) throw new IllegalArgumentException("empty spectrum");

        int zeros = 0;
        for (double x : laplacianEvalsAscending) {
            if (!Double.isFinite(x)) throw new IllegalArgumentException("non-finite eigenvalue: " + x);
            if (Math.abs(x) <= tol) zeros++;
        }
        if (zeros != 1) {
            throw new IllegalArgumentException("Connected Laplacian should have exactly one ~0 eigenvalue; found " + zeros + " (tol=" + tol + ")");
        }
    }

    public static double[] solveSymmetric(DMatrixRMaj A, DMatrixRMaj eigenvectorsColsOut) {
        Objects.requireNonNull(A, "A is null");
        Objects.requireNonNull(eigenvectorsColsOut, "eigenvectorsColsOut is null");

        final int n = A.numRows;
        if (A.numCols != n) throw new IllegalArgumentException("Matrix must be square. Got " + A.numRows + "x" + A.numCols);
        if (n == 0) {
            eigenvectorsColsOut.reshape(0, 0);
            return new double[0];
        }

        // Exactness-preserving validation: do NOT modify A; just reject invalid inputs.
        requireFinite(A);
        requireExactlySymmetric(A);

        // true => compute eigenvectors
        EigenDecomposition_F64<DMatrixRMaj> eig = DecompositionFactory_DDRM.eig(n, true);
        if (!eig.decompose(A.copy())) {
            throw new IllegalStateException("Eigendecomposition failed.");
        }

        double[] vals = new double[n];
        DMatrixRMaj vecs = new DMatrixRMaj(n, n);

        for (int i = 0; i < n; i++) {
            Complex_F64 ev = eig.getEigenvalue(i);
            if (ev == null) throw new IllegalStateException("Missing eigenvalue at i=" + i);

            // For real symmetric matrices, eigenvalues must be real.
            if (ev.getImaginary() != 0.0) {
                throw new IllegalArgumentException(
                        "Matrix is not behaving as real-symmetric: eigenvalue has imaginary part at i=" + i + " imag=" + ev.getImaginary()
                );
            }

            double lambda = ev.getReal();
            if (!Double.isFinite(lambda)) {
                throw new IllegalStateException("Non-finite eigenvalue at i=" + i + ": " + lambda);
            }
            vals[i] = lambda;

            DMatrixRMaj v = eig.getEigenVector(i);
            if (v == null) {
                throw new IllegalStateException("Missing eigenvector at i=" + i);
            }
            if (v.numRows != n || v.numCols != 1) {
                throw new IllegalStateException("Unexpected eigenvector shape at i=" + i + ": " + v.numRows + "x" + v.numCols);
            }

            // Copy v (n x 1) into column i of vecs
            for (int r = 0; r < n; r++) {
                double x = v.get(r, 0);
                if (!Double.isFinite(x)) {
                    throw new IllegalStateException("Non-finite eigenvector entry at i=" + i + " r=" + r + ": " + x);
                }
                vecs.set(r, i, x);
            }
        }

        // Sort eigenpairs by ascending eigenvalue (do NOT rely on EJML order)
        int[] order = argsortDeterministic(vals);

        double[] valsSorted = new double[n];
        eigenvectorsColsOut.reshape(n, n);

        for (int j = 0; j < n; j++) {
            int i = order[j];
            valsSorted[j] = vals[i];

            for (int r = 0; r < n; r++) {
                eigenvectorsColsOut.set(r, j, vecs.get(r, i));
            }

            // Deterministic sign convention:
            // Find entry with largest |value|; force it to be >= 0 by flipping the whole vector if needed.
            normalizeEigenvectorColumnSignInPlace(eigenvectorsColsOut, j);
        }

        return valsSorted;
    }

    public static double[] solveSymmetricEigenvaluesOnly(DMatrixRMaj A) {
        Objects.requireNonNull(A, "A is null");

        final int n = A.numRows;
        if (A.numCols != n) throw new IllegalArgumentException("Matrix must be square. Got " + A.numRows + "x" + A.numCols);
        if (n == 0) return new double[0];

        requireFinite(A);
        requireExactlySymmetric(A);

        EigenDecomposition_F64<DMatrixRMaj> eig = DecompositionFactory_DDRM.eig(n, false);
        if (!eig.decompose(A.copy())) {
            throw new IllegalStateException("Eigendecomposition failed.");
        }

        double[] vals = new double[n];
        for (int i = 0; i < n; i++) {
            Complex_F64 ev = eig.getEigenvalue(i);
            if (ev == null) throw new IllegalStateException("Missing eigenvalue at i=" + i);
            if (ev.getImaginary() != 0.0) {
                throw new IllegalArgumentException(
                        "Matrix is not behaving as real-symmetric: eigenvalue has imaginary part at i=" + i + " imag=" + ev.getImaginary()
                );
            }
            double lambda = ev.getReal();
            if (!Double.isFinite(lambda)) throw new IllegalStateException("Non-finite eigenvalue at i=" + i + ": " + lambda);
            vals[i] = lambda;
        }

        int[] order = argsortDeterministic(vals);
        double[] sorted = new double[n];
        for (int j = 0; j < n; j++) sorted[j] = vals[order[j]];
        return sorted;
    }

    private static void requireFinite(DMatrixRMaj A) {
        final int rows = A.numRows, cols = A.numCols;
        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                double x = A.get(r, c);
                if (!Double.isFinite(x)) {
                    throw new IllegalArgumentException("Matrix contains non-finite entry at (" + r + "," + c + "): " + x);
                }
            }
        }
    }

    private static void requireExactlySymmetric(DMatrixRMaj A) {
        final int n = A.numRows;
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                double aij = A.get(i, j);
                double aji = A.get(j, i);
                if (Double.doubleToLongBits(aij) != Double.doubleToLongBits(aji)) {
                    throw new IllegalArgumentException(
                            "Matrix must be exactly symmetric: A(" + i + "," + j + ")=" + aij +
                                    " != A(" + j + "," + i + ")=" + aji
                    );
                }
            }
        }
    }

    private static void normalizeEigenvectorColumnSignInPlace(DMatrixRMaj V, int col) {
        final int n = V.numRows;

        int argMax = 0;
        double maxAbs = 0.0;

        for (int r = 0; r < n; r++) {
            double x = V.get(r, col);
            double ax = Math.abs(x);
            if (ax > maxAbs) {
                maxAbs = ax;
                argMax = r;
            }
        }

        // If the vector is all zeros (shouldn't happen), do nothing.
        if (maxAbs == 0.0) return;

        double pivot = V.get(argMax, col);
        if (pivot < 0.0) {
            for (int r = 0; r < n; r++) {
                V.set(r, col, -V.get(r, col));
            }
        }
    }

    private static int[] argsortDeterministic(double[] a) {
        final int n = a.length;
        int[] idx = new int[n];
        for (int i = 0; i < n; i++) idx[i] = i;
        quicksortIdx(idx, a, 0, n - 1);
        return idx;
    }

    private static void quicksortIdx(int[] idx, double[] a, int lo, int hi) {
        while (lo < hi) {
            int p = partition(idx, a, lo, hi);
            if (p - lo < hi - p) {
                quicksortIdx(idx, a, lo, p - 1);
                lo = p + 1;
            } else {
                quicksortIdx(idx, a, p + 1, hi);
                hi = p - 1;
            }
        }
    }

    private static int partition(int[] idx, double[] a, int lo, int hi) {
        int mid = lo + ((hi - lo) >>> 1);
        int pivotIndex = medianOf3(idx, a, lo, mid, hi);
        swap(idx, pivotIndex, hi);

        int pivotId = idx[hi];
        double pivotVal = a[pivotId];

        int i = lo;
        for (int j = lo; j < hi; j++) {
            int idJ = idx[j];
            if (lessByValueThenIndex(a[idJ], idJ, pivotVal, pivotId)) {
                swap(idx, i, j);
                i++;
            }
        }
        swap(idx, i, hi);
        return i;
    }

    private static int medianOf3(int[] idx, double[] a, int i, int j, int k) {
        int ai = idx[i], aj = idx[j], ak = idx[k];
        double vi = a[ai], vj = a[aj], vk = a[ak];

        if (lessByValueThenIndex(vi, ai, vj, aj)) {
            if (lessByValueThenIndex(vj, aj, vk, ak)) return j;
            if (lessByValueThenIndex(vi, ai, vk, ak)) return k;
            return i;
        } else {
            if (lessByValueThenIndex(vi, ai, vk, ak)) return i;
            if (lessByValueThenIndex(vj, aj, vk, ak)) return k;
            return j;
        }
    }

    private static boolean lessByValueThenIndex(double v1, int i1, double v2, int i2) {
        int c = Double.compare(v1, v2);
        if (c != 0) return c < 0;
        return i1 < i2;
    }

    private static void swap(int[] a, int i, int j) {
        int t = a[i]; a[i] = a[j]; a[j] = t;
    }
}
