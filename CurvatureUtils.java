package iso;

import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.factory.LinearSolverFactory_DDRM;
import org.ejml.interfaces.linsol.LinearSolverDense;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public final class CurvatureUtils {

    private CurvatureUtils() {}

    public enum SamplingDomain {
        LOG_T,
        LINEAR_T
    }

    private enum TimeGridMode {
        GEOMSPACE,
        CHEBYSHEV_NODES_LOGT
    }

    public static final class Result {
        /** Time grid t[j], length T. */
        public final double[] t;

        /** Raw heat-kernel diagonal K[u][j], shape (n,T). */
        public final double[][] Kdiag;

        /** Coefficients per vertex (shape depends on includeC0 / degrees). */
        public final double[][] coeffs;

        /** Quantized coefficients, same shape as coeffs. */
        public final long[][] qCoeffs;

        /** Quantization scale used. */
        public final double coeffQuantScale;

        /** Per-vertex GOF computed in the fitted space (RAW-K or LOG-scaled, depending on Constants). */
        public final double[] r2;
        public final double[] adjR2;
        public final double[] rmse;
        public final double[] sse;

        private Result(
                double[] t,
                double[][] Kdiag,
                double[][] coeffs,
                long[][] qCoeffs,
                double coeffQuantScale,
                double[] r2,
                double[] adjR2,
                double[] rmse,
                double[] sse
        ) {
            this.t = t;
            this.Kdiag = Kdiag;
            this.coeffs = coeffs;
            this.qCoeffs = qCoeffs;
            this.coeffQuantScale = coeffQuantScale;
            this.r2 = r2;
            this.adjR2 = adjR2;
            this.rmse = rmse;
            this.sse = sse;
        }
    }

    public static void printRawCoeffs(String label, double[][] coeffs, int... vertices)
    {
    	System.out.println("\nRAW COEFFS: " + label);

    	for (int v : vertices) {
    	System.out.println("v=" + v + " -> " +
    	Arrays.toString(coeffs[v]));
    	}
    }

    /**
     * End-to-end:
     *  1) auto choose T based on max degree
     *  2) build time grid from effective spectrum mu = D*lambda^s
     *  3) compute rotation-invariant heat-kernel diagonal
     *  4) fit coefficients (RAW or LOG space controlled by Constants)
     *  5) compute per-vertex GOF in fitted space
     *  6) quantize coefficients
     */
    public static Result compute(
            double[] laplacianEigenvaluesAscending,   // (n,)
            DMatrixRMaj eigenvectorsCols,            // (n,n), columns aligned with ascending eigenvalues
            double ds,
            int[] degreesOrNull,
            int defaultDegree,
            boolean includeC0,
            boolean clipNonFiniteKToZero,
            SamplingDomain samplingDomain,
            double s,
            double D,
            double coeffQuantScale
    ) {
        if (samplingDomain == null) samplingDomain = SamplingDomain.LOG_T;

        int n = laplacianEigenvaluesAscending.length;
        if (eigenvectorsCols.numRows != n || eigenvectorsCols.numCols != n) {
            throw new IllegalArgumentException("eigenvectorsCols must be n x n");
        }
        if (!Double.isFinite(ds)) throw new IllegalArgumentException("ds must be finite");
        if (!Double.isFinite(s) || s <= 0.0) throw new IllegalArgumentException("s must be finite and > 0");
        if (!Double.isFinite(D) || D <= 0.0) throw new IllegalArgumentException("D must be finite and > 0");

        int maxDeg = maxDegree(degreesOrNull, n, defaultDegree);
        int pMax = includeC0 ? (maxDeg + 1) : maxDeg;
        int T = chooseTimePointsFromP(pMax);

        // Infer tMin/tMax from effective spectrum
        double[] tMinMax = inferTMinMaxFromEigenvalues(
                laplacianEigenvaluesAscending,
                Constants.CURV_TGRID_EIG_TOL,
                Constants.CURV_TGRID_SPAN_MAX,
                Constants.CURV_TGRID_MAX_LAM_T,
                Constants.CURV_TGRID_MIN_RATIO,
                s,
                D
        );
        
        //tMinMax[1] = 85.0;
        //System.out.println("tmin " + tMinMax[0]);
        //System.out.println("tmax " + tMinMax[1]);
        
        TimeGridMode gridMode =
                (Constants.CURV_TGRID_USE_CHEBYSHEV_NODES_LOGT && samplingDomain == SamplingDomain.LOG_T)
                        ? TimeGridMode.CHEBYSHEV_NODES_LOGT
                        : TimeGridMode.GEOMSPACE;

        double[] t = (gridMode == TimeGridMode.CHEBYSHEV_NODES_LOGT)
                ? chebyshevNodesTimeGridLogT(tMinMax[0], tMinMax[1], T)
                : geomspace(tMinMax[0], tMinMax[1], T);

        double[][] Kdiag = heatKernelDiagGrouped(
                laplacianEigenvaluesAscending,
                eigenvectorsCols,
                t,
                Constants.CURV_MU_GROUP_ABS_TOL,
                Constants.CURV_MU_GROUP_REL_TOL,
                Constants.CURV_MU_GROUP_ULP_MULT,
                s,
                D
        );

        FitResult fit = Constants.CURV_FIT_IN_LOG_SPACE
                ? fit_LOGSCALED(Kdiag, t, ds, degreesOrNull, defaultDegree, includeC0, clipNonFiniteKToZero, samplingDomain)
                : fit_RAW(Kdiag, t, ds, degreesOrNull, defaultDegree, includeC0, clipNonFiniteKToZero, samplingDomain);

        long[][] q = quantize2D(fit.coeffs, coeffQuantScale);

        return new Result(
                t,
                Kdiag,
                fit.coeffs,
                q,
                coeffQuantScale,
                fit.r2,
                fit.adjR2,
                fit.rmse,
                fit.sse
        );
    }

    /** Convenience overload using project defaults for s, D, and coeff quant scale. */
    public static Result compute(
            double[] laplacianEigenvaluesAscending,
            DMatrixRMaj eigenvectorsCols,
            double ds,
            int defaultDegree
    ) {
        return compute(
                laplacianEigenvaluesAscending,
                eigenvectorsCols,
                ds,
                null,
                defaultDegree,
                true,
                true,
                SamplingDomain.LOG_T,
                Constants.SPECTRAL_EXPONENT,
                Constants.DIFFUSION_CONSTANT,
                Constants.CURV_COEFF_DEFAULT_SCALE
        );
    }

    private static double[] inferTMinMaxFromEigenvalues(
            double[] eigenvalues,
            double eigTol,
            double tSpanMax,
            double maxLamT,
            double minRatio,
            double s,
            double D
    ) {
        if (eigenvalues == null || eigenvalues.length == 0)
            throw new IllegalArgumentException("eigenvalues empty");

        // Treat tiny eigenvalues as numerical zero deterministically.
        // This fixes permutation-dependent solver noise like +8e-16 vs -2e-16 for the zero mode.
        final double zeroLam = Constants.EIGEN_ZERO_CLAMP;

        double muMin = Double.POSITIVE_INFINITY;
        double muMax = 0.0;
        boolean any = false;

        for (double lamRaw : eigenvalues) {
            if (!Double.isFinite(lamRaw)) continue;

            double lam = lamRaw;

            // NEW: hard clamp around zero (both signs)
            if (Math.abs(lam) <= zeroLam) lam = 0.0;

            // Keep original safety
            if (lam < 0.0) lam = 0.0;

            // Ignore true/near-zero eigenvalues
            if (lam == 0.0) continue;

            double mu = effectiveMu(lam, s, D);

            // Also protect against tiny mu noise (eigTol is in mu-space)
            if (mu > eigTol) {
                any = true;
                if (mu < muMin) muMin = mu;
                if (mu > muMax) muMax = mu;
            }
        }

        if (!any)
            throw new IllegalArgumentException("No positive eigenvalues above eigTol; cannot build time grid.");

        double tMin = 1.0 / muMax;
        double tMax = 1.0 / muMin;

        tMax = Math.min(tMax, tMin * tSpanMax);
        tMax = Math.min(tMax, maxLamT / muMax);

        //System.out.println("muMax: " + muMax);
        //System.out.println("muMin: " + muMin);
        
        if (tMax <= tMin * minRatio) tMax = tMin * minRatio;

        if (!(Double.isFinite(tMin) && Double.isFinite(tMax)) || tMin <= 0.0 || tMax <= 0.0) {
            throw new IllegalArgumentException("Invalid tMin/tMax");
        }
        if (tMax <= tMin) tMax = tMin * minRatio;

        return new double[]{tMin, tMax};
    }

    private static double[] geomspace(double a, double b, int n) {
        if (n < 2) throw new IllegalArgumentException("n must be >= 2");
        double[] out = new double[n];
        double logA = Math.log(a);
        double logB = Math.log(b);
        for (int i = 0; i < n; i++) {
            double tt = (double) i / (double) (n - 1);
            out[i] = Math.exp(logA + tt * (logB - logA));
        }
        for (int i = 1; i < n; i++) {
            if (!(out[i] > out[i - 1])) out[i] = Math.nextUp(out[i - 1]);
        }
        return out;
    }

    /** Next-trick time grid: Chebyshev nodes in log(t), mapped back to t. */
    private static double[] chebyshevNodesTimeGridLogT(double tMin, double tMax, int nT) {
        if (nT < 2) throw new IllegalArgumentException("nT must be >= 2");
        if (!(tMin > 0.0 && tMax > 0.0 && tMax > tMin)) throw new IllegalArgumentException("invalid tMin/tMax");

        double lo = Math.log(tMin);
        double hi = Math.log(tMax);

        double[] t = new double[nT];
        for (int j = 0; j < nT; j++) {
            double x = Math.cos(Math.PI * (j + 0.5) / nT); // in (1..-1)
            double logt = 0.5 * (hi + lo) + 0.5 * (hi - lo) * x;
            t[j] = Math.exp(logt);
        }
        Arrays.sort(t);
        for (int i = 1; i < nT; i++) {
            if (!(t[i] > t[i - 1])) t[i] = Math.nextUp(t[i - 1]);
        }
        return t;
    }

    public static double[][] heatKernelDiagGrouped(
            double[] eigenvaluesAscending,
            DMatrixRMaj eigenvectorsCols,
            double[] t,
            double absTol,
            double relTol,
            double ulpMult,
            double s,
            double D
    ) {
        int n = eigenvaluesAscending.length;
        if (n == 0) throw new IllegalArgumentException("empty spectrum");
        if (eigenvectorsCols.numRows != n || eigenvectorsCols.numCols != n) {
            throw new IllegalArgumentException("eigenvectorsCols must be n x n");
        }
        if (t == null) throw new IllegalArgumentException("t is null");
        int T = t.length;
        if (T == 0) return new double[n][0];

        for (double tj : t) {
            if (!Double.isFinite(tj) || tj <= 0.0) throw new IllegalArgumentException("t must be finite and > 0");
        }
        if (!(absTol >= 0.0) || !(relTol >= 0.0) || !(ulpMult >= 0.0)) {
            throw new IllegalArgumentException("absTol/relTol/ulpMult must be >= 0");
        }

        // mu[i] = D * max(lambda,0)^s
        double[] mu = new double[n];
        for (int i = 0; i < n; i++) {
            double lam = eigenvaluesAscending[i];
            if (!Double.isFinite(lam)) {
                mu[i] = lam;
                continue;
            }
            if (lam < 0.0) lam = 0.0;
            mu[i] = (lam == 0.0) ? 0.0 : effectiveMu(lam, s, D);
        }

        // Group consecutive near-equal mu (local scale + ulp floor)
        List<int[]> groups = new ArrayList<>();
        int a = 0;
        for (int i = 1; i < n; i++) {
            double m0 = mu[i - 1];
            double m1 = mu[i];
            double diff = Math.abs(m1 - m0);

            double scale = Math.max(1.0, Math.max(Math.abs(m0), Math.abs(m1)));
            double tol = Math.max(absTol, relTol * scale);
            tol = Math.max(tol, ulpMult * Math.ulp(scale));

            if (!(diff <= tol)) {
                groups.add(new int[]{a, i});
                a = i;
            }
        }
        groups.add(new int[]{a, n});

        // Precompute projection masses per group: P_g(u) = sum_{k in group} v_{u,k}^2
        int gCount = groups.size();
        double[][] proj = new double[gCount][n];
        double[] muRep = new double[gCount];

        for (int gi = 0; gi < gCount; gi++) {
            int start = groups.get(gi)[0];
            int end = groups.get(gi)[1];

            double sum = 0.0;
            int cnt = 0;
            for (int k = start; k < end; k++) {
                double m = mu[k];
                if (Double.isFinite(m)) { sum += m; cnt++; }
            }
            muRep[gi] = (cnt == 0) ? Double.NaN : (sum / (double) cnt);

            double[] pg = proj[gi];
            for (int u = 0; u < n; u++) {
                double ss = 0.0;
                for (int k = start; k < end; k++) {
                    double v = eigenvectorsCols.get(u, k);
                    ss += v * v;
                }
                pg[u] = ss;
            }
        }

        // Assemble K(u,t)
        double[][] K = new double[n][T];
        for (int j = 0; j < T; j++) {
            double tj = t[j];

            for (int gi = 0; gi < gCount; gi++) {
                double m = muRep[gi];
                if (!Double.isFinite(m)) continue;

                double arg = -tj * m;
                if (arg < -Constants.CURV_EXP_CLIP) continue; // underflow -> 0

                double w = Math.exp(arg);
                double[] pg = proj[gi];
                for (int u = 0; u < n; u++) {
                    K[u][j] += w * pg[u];
                }
            }
        }
        return K;
    }

    private static final class FitResult {
        final double[][] coeffs;
        final double[] r2;
        final double[] adjR2;
        final double[] rmse;
        final double[] sse;

        FitResult(double[][] coeffs, double[] r2, double[] adjR2, double[] rmse, double[] sse) {
            this.coeffs = coeffs;
            this.r2 = r2;
            this.adjR2 = adjR2;
            this.rmse = rmse;
            this.sse = sse;
        }
    }

    /**
     * RAW fit:
     *   K(u,t) ≈ pref(t) * Σ a_k T_k(x(t)), where pref(t)=(4πt)^(-ds/2)
     * GOF computed in K-space.
     */
    private static FitResult fit_RAW(
            double[][] K,
            double[] t,
            double ds,
            int[] degreesOrNull,
            int defaultDegree,
            boolean includeC0,
            boolean clipNonFiniteToZero,
            SamplingDomain samplingDomain
    ) {
        RectK rk = sanitizeK(K, t, clipNonFiniteToZero);
        double[][] Kuse = rk.K;
        int n = rk.n;
        int T = rk.T;

        int[] degs = makeDegrees(degreesOrNull, n, defaultDegree);
        int maxDeg = maxIn(degs);

        double[] x = mapToMinus1Plus1(t, samplingDomain);

        // pref(t) = (4πt)^(-ds/2) in log-space with clipping
        double[] pref = new double[T];
        for (int j = 0; j < T; j++) {
            double logPref = (-ds / 2.0) * Math.log(4.0 * Math.PI * t[j]);
            if (logPref < -Constants.CURV_LOG_EXP_CLIP) logPref = -Constants.CURV_LOG_EXP_CLIP;
            if (logPref >  Constants.CURV_LOG_EXP_CLIP) logPref =  Constants.CURV_LOG_EXP_CLIP;
            pref[j] = Math.exp(logPref);
        }

        DMatrixRMaj Phi = chebyshevTMatrix(x, maxDeg);

        // A = pref * Phi  => (T, maxDeg+1)
        DMatrixRMaj A = new DMatrixRMaj(T, maxDeg + 1);
        for (int r = 0; r < T; r++) {
            double pr = pref[r];
            for (int c = 0; c <= maxDeg; c++) {
                A.set(r, c, Phi.get(r, c) * pr);
            }
        }

        // Y = K^T  => (T,n)
        DMatrixRMaj Y = new DMatrixRMaj(T, n);
        for (int u = 0; u < n; u++) {
            for (int j = 0; j < T; j++) {
                Y.set(j, u, Kuse[u][j]);
            }
        }

        // B = pinv(A) * Y   => (maxDeg+1, n)
        //DMatrixRMaj pinvA = new DMatrixRMaj(maxDeg + 1, T);
        //CommonOps_DDRM.pinv(A, pinvA);

        //DMatrixRMaj B = new DMatrixRMaj(maxDeg + 1, n);
        //CommonOps_DDRM.mult(pinvA, Y, B);
        
        // Solve least-squares: A * B = Y   (A is T x (maxDeg+1), Y is T x n)
        DMatrixRMaj B = solveLeastSquaresQR(A, Y);
        
        double[][] coeffs = unpackCoeffs(B, degs, maxDeg, includeC0);

        // GOF in K space
        double[] r2 = new double[n];
        double[] adjR2 = new double[n];
        double[] rmse = new double[n];
        double[] sse = new double[n];

        for (int u = 0; u < n; u++) {
            int d = degs[u];
            int p = includeC0 ? (d + 1) : d;

            // mean of K(u,:)
            double mean = 0.0;
            for (int j = 0; j < T; j++) mean += Kuse[u][j];
            mean /= (double) T;

            double ssTot = 0.0;
            double ssRes = 0.0;

            for (int j = 0; j < T; j++) {
                double yObs = Kuse[u][j];
                double yHat = predictK_RAW(u, j, d, includeC0, coeffs, Phi, pref);

                double r = yObs - yHat;
                ssRes += r * r;

                double dy = yObs - mean;
                ssTot += dy * dy;
            }

            sse[u] = ssRes;
            rmse[u] = Math.sqrt(ssRes / (double) T);

            double r2u;
            if (ssTot == 0.0) {
                r2u = (ssRes == 0.0) ? 1.0 : 0.0;
            } else {
                r2u = 1.0 - (ssRes / ssTot);
                if (r2u < -1e-12) r2u = -1e-12;
                if (r2u > 1.0 + 1e-12) r2u = 1.0;
            }
            r2[u] = r2u;

            if (T <= p + 1) adjR2[u] = Double.NaN;
            else adjR2[u] = 1.0 - (1.0 - r2u) * (T - 1.0) / (T - p - 1.0);
        }

        return new FitResult(coeffs, r2, adjR2, rmse, sse);
    }

    /**
     * LOG-scaled fit (recommended):
     *   y(u,t) = log(max(K(u,t), eps)) + (ds/2)*log(4πt)
     *   y ≈ Σ a_k T_k(x(t))
     * GOF computed in y space.
     */
    private static FitResult fit_LOGSCALED(
            double[][] K,
            double[] t,
            double ds,
            int[] degreesOrNull,
            int defaultDegree,
            boolean includeC0,
            boolean clipNonFiniteToZero,
            SamplingDomain samplingDomain
    ) {
        RectK rk = sanitizeK(K, t, clipNonFiniteToZero);
        double[][] Kuse = rk.K;
        int n = rk.n;
        int T = rk.T;

        int[] degs = makeDegrees(degreesOrNull, n, defaultDegree);
        int maxDeg = maxIn(degs);

        double[] x = mapToMinus1Plus1(t, samplingDomain);
        DMatrixRMaj Phi = chebyshevTMatrix(x, maxDeg);

        // Build Y = y^T (T,n) where y = log(K) + (ds/2)log(4πt)
        DMatrixRMaj Y = new DMatrixRMaj(T, n);
        double halfDs = ds / 2.0;

        for (int u = 0; u < n; u++) {
            for (int j = 0; j < T; j++) {
                double kij = Kuse[u][j];
                if (!Double.isFinite(kij) || kij <= 0.0) kij = Constants.CURV_LOG_FIT_MIN_K;

                double y = Math.log(kij) + halfDs * Math.log(4.0 * Math.PI * t[j]);

                // optional clip
                if (y < -Constants.CURV_LOG_EXP_CLIP) y = -Constants.CURV_LOG_EXP_CLIP;
                if (y >  Constants.CURV_LOG_EXP_CLIP) y =  Constants.CURV_LOG_EXP_CLIP;

                Y.set(j, u, y);
            }
        }

        // B = pinv(Phi) * Y  => (maxDeg+1, n)
        //DMatrixRMaj pinvPhi = new DMatrixRMaj(maxDeg + 1, T);
        //CommonOps_DDRM.pinv(Phi, pinvPhi);

        //DMatrixRMaj B = new DMatrixRMaj(maxDeg + 1, n);
        //CommonOps_DDRM.mult(pinvPhi, Y, B);

        // Solve least-squares: Phi * B = Y  (Phi is T x (maxDeg+1), Y is T x n)
        DMatrixRMaj B = solveLeastSquaresQR(Phi, Y);

        double[][] coeffs = unpackCoeffs(B, degs, maxDeg, includeC0);

        // GOF in y-space
        double[] r2 = new double[n];
        double[] adjR2 = new double[n];
        double[] rmse = new double[n];
        double[] sse = new double[n];

        for (int u = 0; u < n; u++) {
            int d = degs[u];
            int p = includeC0 ? (d + 1) : d;

            double mean = 0.0;
            for (int j = 0; j < T; j++) mean += Y.get(j, u);
            mean /= (double) T;

            double ssTot = 0.0;
            double ssRes = 0.0;

            for (int j = 0; j < T; j++) {
                double yObs = Y.get(j, u);
                double yHat = predictY_LOG(u, j, d, includeC0, coeffs, Phi);

                double r = yObs - yHat;
                ssRes += r * r;

                double dy = yObs - mean;
                ssTot += dy * dy;
            }

            sse[u] = ssRes;
            rmse[u] = Math.sqrt(ssRes / (double) T);

            double r2u;
            if (ssTot == 0.0) {
                r2u = (ssRes == 0.0) ? 1.0 : 0.0;
            } else {
                r2u = 1.0 - (ssRes / ssTot);
                if (r2u < -1e-12) r2u = -1e-12;
                if (r2u > 1.0 + 1e-12) r2u = 1.0;
            }
            r2[u] = r2u;

            if (T <= p + 1) adjR2[u] = Double.NaN;
            else adjR2[u] = 1.0 - (1.0 - r2u) * (T - 1.0) / (T - p - 1.0);
        }

        return new FitResult(coeffs, r2, adjR2, rmse, sse);
    }
    
    private static double predictK_RAW(
            int u, int j, int degreeU, boolean includeC0,
            double[][] coeffs, DMatrixRMaj Phi, double[] pref
    ) {
        double scale = pref[j];
        double sum = 0.0;
        if (includeC0) {
            for (int k = 0; k <= degreeU; k++) sum += coeffs[u][k] * (Phi.get(j, k) * scale);
        } else {
            for (int k = 1; k <= degreeU; k++) sum += coeffs[u][k - 1] * (Phi.get(j, k) * scale);
        }
        return sum;
    }

    private static double predictY_LOG(
            int u, int j, int degreeU, boolean includeC0,
            double[][] coeffs, DMatrixRMaj Phi
    ) {
        double sum = 0.0;
        if (includeC0) {
            for (int k = 0; k <= degreeU; k++) sum += coeffs[u][k] * Phi.get(j, k);
        } else {
            for (int k = 1; k <= degreeU; k++) sum += coeffs[u][k - 1] * Phi.get(j, k);
        }
        return sum;
    }

    /**
     * Solve least-squares min ||A*X - B|| using QR (stable).
     * A is (m x n), B is (m x k), returns X as (n x k).
     */
    private static DMatrixRMaj solveLeastSquaresQR(DMatrixRMaj A, DMatrixRMaj B) {
        if (A.numRows != B.numRows) {
            throw new IllegalArgumentException("Row mismatch: A is " + A.numRows + "x" + A.numCols +
                    " but B is " + B.numRows + "x" + B.numCols);
        }

        // QR least squares solver
        LinearSolverDense<DMatrixRMaj> solver = LinearSolverFactory_DDRM.leastSquares(A.numRows, A.numCols);

        if (!solver.setA(A)) {
            throw new IllegalStateException("QR solver failed to setA (A may be singular or ill-conditioned).");
        }

        DMatrixRMaj X = new DMatrixRMaj(A.numCols, B.numCols);
        solver.solve(B, X);
        return X;
    }

    private static final class RectK {
        final double[][] K;
        final int n;
        final int T;
        RectK(double[][] K, int n, int T) { this.K = K; this.n = n; this.T = T; }
    }

    private static RectK sanitizeK(double[][] K, double[] t, boolean clipNonFiniteToZero) {
        if (K == null || t == null) throw new IllegalArgumentException("K/t null");
        int n = K.length;
        int T = (n == 0) ? 0 : K[0].length;
        if (t.length != T) throw new IllegalArgumentException("t length must equal K columns (T)");
        if (T == 0) return new RectK(new double[n][0], n, 0);

        // validate t (strictly increasing)
        for (int i = 0; i < T; i++) {
            double ti = t[i];
            if (!Double.isFinite(ti) || ti <= 0.0) throw new IllegalArgumentException("t must be finite and > 0");
            if (i > 0 && !(t[i] > t[i - 1])) throw new IllegalArgumentException("t should be strictly increasing");
        }

        if (clipNonFiniteToZero) {
            double[][] out = new double[n][T];
            for (int u = 0; u < n; u++) {
                if (K[u].length != T) throw new IllegalArgumentException("K must be rectangular (n,T)");
                for (int j = 0; j < T; j++) {
                    double v = K[u][j];
                    out[u][j] = Double.isFinite(v) ? v : 0.0;
                }
            }
            return new RectK(out, n, T);
        } else {
            for (int u = 0; u < n; u++) {
                if (K[u].length != T) throw new IllegalArgumentException("K must be rectangular (n,T)");
                for (int j = 0; j < T; j++) {
                    if (!Double.isFinite(K[u][j])) {
                        throw new IllegalArgumentException("K contains non-finite; enable clipNonFiniteKToZero");
                    }
                }
            }
            return new RectK(K, n, T);
        }
    }

    private static int[] makeDegrees(int[] degreesOrNull, int n, int defaultDegree) {
        int[] degs = new int[n];
        if (degreesOrNull == null) {
            if (defaultDegree < 0) throw new IllegalArgumentException("defaultDegree must be >= 0");
            Arrays.fill(degs, defaultDegree);
            return degs;
        }
        if (degreesOrNull.length != n) throw new IllegalArgumentException("degrees length mismatch");
        for (int i = 0; i < n; i++) {
            int d = degreesOrNull[i];
            if (d < 0) throw new IllegalArgumentException("degrees must be >= 0");
            degs[i] = d;
        }
        return degs;
    }

    private static int maxIn(int[] a) {
        int m = 0;
        for (int v : a) m = Math.max(m, v);
        return m;
    }

    private static int maxDegree(int[] degreesOrNull, int n, int defaultDegree) {
        if (n == 0) return 0;
        if (degreesOrNull == null) return defaultDegree;
        if (degreesOrNull.length != n) throw new IllegalArgumentException("degrees length mismatch");
        int m = 0;
        for (int d : degreesOrNull) m = Math.max(m, d);
        return m;
    }

    /** Choose T based on number of parameters p. */
    public static int chooseTimePointsFromP(int p) {
        if (p < 1) p = 1;
        int T = Constants.CURV_TGRID_POINTS_PER_COEFF * p;
        if (T < Constants.CURV_TGRID_MIN_POINTS) T = Constants.CURV_TGRID_MIN_POINTS;
        if (T > Constants.CURV_TGRID_MAX_POINTS) T = Constants.CURV_TGRID_MAX_POINTS;
        return T;
    }

    private static double[] mapToMinus1Plus1(double[] t, SamplingDomain domain) {
        int T = t.length;
        double[] z = new double[T];

        if (domain == SamplingDomain.LOG_T) {
            double lo = Double.POSITIVE_INFINITY, hi = Double.NEGATIVE_INFINITY;
            for (int j = 0; j < T; j++) {
                double v = Math.log(t[j]);
                z[j] = v;
                lo = Math.min(lo, v);
                hi = Math.max(hi, v);
            }
            return affineToMinus1Plus1(z, lo, hi);
        }

        if (domain == SamplingDomain.LINEAR_T) {
            double lo = Double.POSITIVE_INFINITY, hi = Double.NEGATIVE_INFINITY;
            for (int j = 0; j < T; j++) {
                double v = t[j];
                z[j] = v;
                lo = Math.min(lo, v);
                hi = Math.max(hi, v);
            }
            return affineToMinus1Plus1(z, lo, hi);
        }

        throw new IllegalArgumentException("Unknown SamplingDomain: " + domain);
    }

    private static double[] affineToMinus1Plus1(double[] z, double lo, double hi) {
        int T = z.length;
        double[] x = new double[T];

        if (hi <= lo + 1e-300) {
            Arrays.fill(x, 0.0);
            return x;
        }
        for (int j = 0; j < T; j++) {
            x[j] = (2.0 * (z[j] - lo) / (hi - lo)) - 1.0;
        }
        return x;
    }

    private static DMatrixRMaj chebyshevTMatrix(double[] x, int maxDeg) {
        int T = x.length;
        DMatrixRMaj Phi = new DMatrixRMaj(T, maxDeg + 1);

        for (int r = 0; r < T; r++) Phi.set(r, 0, 1.0);
        if (maxDeg == 0) return Phi;

        for (int r = 0; r < T; r++) Phi.set(r, 1, x[r]);

        for (int k = 2; k <= maxDeg; k++) {
            for (int r = 0; r < T; r++) {
                double Tk1 = Phi.get(r, k - 1);
                double Tk2 = Phi.get(r, k - 2);
                Phi.set(r, k, 2.0 * x[r] * Tk1 - Tk2);
            }
        }
        return Phi;
    }

    private static double[][] unpackCoeffs(DMatrixRMaj B, int[] degs, int maxDeg, boolean includeC0) {
        int n = degs.length;

        if (includeC0) {
            double[][] out = new double[n][maxDeg + 1];
            for (int u = 0; u < n; u++) {
                int d = degs[u];
                for (int k = 0; k <= d; k++) out[u][k] = B.get(k, u);
            }
            return out;
        } else {
            double[][] out = new double[n][maxDeg];
            for (int u = 0; u < n; u++) {
                int d = degs[u];
                for (int k = 1; k <= d; k++) out[u][k - 1] = B.get(k, u);
            }
            return out;
        }
    }

    public static long[][] quantize2D(double[][] a, double scale) {
        if (scale <= 0) throw new IllegalArgumentException("scale must be > 0");
        int n = a.length;
        long[][] q = new long[n][];

        final double clamp = 0.5 / scale;
        final long B = (long) Constants.CURV_COEFF_QUANT_BUCKET; // make sure this exists and is long
        final long half = B / 2;

        for (int i = 0; i < n; i++) {
            double[] row = a[i];
            q[i] = new long[row.length];
            for (int j = 0; j < row.length; j++) {
                double x = row[j];
                if (!Double.isFinite(x)) {
                    q[i][j] = 0L;
                    continue;
                }
                if (Math.abs(x) < clamp) x = 0.0;

                long raw = Math.round(x * scale);

                if (B > 1) {
                    // symmetric snapping to nearest multiple of B (works for negative too)
                    raw = (raw >= 0)
                            ? ((raw + half) / B) * B
                            : ((raw - half) / B) * B;
                }

                q[i][j] = raw;
            }
        }
        return q;
    }

    private static double effectiveMu(double lambdaPositive, double s, double D) {
        if (!(lambdaPositive > 0.0) || !Double.isFinite(lambdaPositive)) return 0.0;

        if (s == 1.0) {
            double mu = D * lambdaPositive;
            return Double.isFinite(mu) ? mu : Double.POSITIVE_INFINITY;
        }

        double logMu = Math.log(D) + s * Math.log(lambdaPositive);
        if (logMu < Constants.CURV_LOG_MU_MIN) logMu = Constants.CURV_LOG_MU_MIN;
        if (logMu > Constants.CURV_LOG_MU_MAX) logMu = Constants.CURV_LOG_MU_MAX;
        return Math.exp(logMu);
    }
    
    public static final class R2Stats {
        public final double min;
        public final double median;
        public final double mean;
        public final double max;
        public final int count;

        public R2Stats(double min, double median, double mean, double max, int count) {
            this.min = min;
            this.median = median;
            this.mean = mean;
            this.max = max;
            this.count = count;
        }

        @Override
        public String toString() {
            return String.format(
                    "R2Stats{min=%.6f, median=%.6f, mean=%.6f, max=%.6f, n=%d}",
                    min, median, mean, max, count
            );
        }
    }

    public static R2Stats computeR2Stats(double[] r2) {
        if (r2 == null || r2.length == 0)
            throw new IllegalArgumentException("r2 empty");

        int n = r2.length;

        // copy so we don't mutate original
        double[] tmp = new double[n];
        int k = 0;

        double min = Double.POSITIVE_INFINITY;
        double max = Double.NEGATIVE_INFINITY;
        double sum = 0.0;

        for (double v : r2) {
            if (!Double.isFinite(v)) continue;

            tmp[k++] = v;
            sum += v;

            if (v < min) min = v;
            if (v > max) max = v;
        }

        if (k == 0)
            throw new IllegalArgumentException("no finite r2 values");

        Arrays.sort(tmp, 0, k);

        double median;
        if ((k & 1) == 0) {
            median = 0.5 * (tmp[k/2 - 1] + tmp[k/2]);
        } else {
            median = tmp[k/2];
        }

        double mean = sum / k;

        return new R2Stats(min, median, mean, max, k);
    }
}
