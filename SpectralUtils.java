package iso;

import java.util.Arrays;
import java.util.Objects;

public final class SpectralUtils {

    private SpectralUtils() {}

    /** Fit-space choice. */
    public enum FitSpace {
        /** x=log(lambda), y=log(rho). Only here Ds=2*slope is meaningful. */
        LOG_LOG,
        /** x=lambda, y=log(rho). Ds not defined. */
        LINEAR_LOG,
        /** x=lambda, y=rho. Ds not defined. */
        LINEAR_LINEAR
    }

    /** Result bundle. */
    public static final class DsResult {
        public final double ds;      // meaningful only for LOG_LOG; else NaN
        public final double gof;     // adjusted R^2 for chosen window
        public final double slope;
        public final double intercept;
        public final FitSpace fitSpace;
        public final int startPointIndex;
        public final int endPointIndexExclusive;
        public final int usedPoints;

        private DsResult(
                double ds,
                double gof,
                double slope,
                double intercept,
                FitSpace fitSpace,
                int startPointIndex,
                int endPointIndexExclusive
        ) {
            this.ds = ds;
            this.gof = gof;
            this.slope = slope;
            this.intercept = intercept;
            this.fitSpace = fitSpace;
            this.startPointIndex = startPointIndex;
            this.endPointIndexExclusive = endPointIndexExclusive;
            this.usedPoints = Math.max(0, endPointIndexExclusive - startPointIndex);
        }

        @Override
        public String toString() {
            return "DsResult{" +
                    "ds=" + ds +
                    ", gof=" + gof +
                    ", slope=" + slope +
                    ", intercept=" + intercept +
                    ", fitSpace=" + fitSpace +
                    ", window=[" + startPointIndex + "," + endPointIndexExclusive + ")" +
                    ", usedPoints=" + usedPoints +
                    '}';
        }
    }

    /** Default: LOG_LOG using Constants, with (s=1, D=1). */
    public static DsResult estimateDs(double[] eigenvalues) {
        return estimateDs(eigenvalues, FitSpace.LOG_LOG, 1.0, 1.0);
    }

    /** Default: LOG_LOG using Constants, with given s (D=1). */
    public static DsResult estimateDs(double[] eigenvalues, double s) {
        return estimateDs(eigenvalues, FitSpace.LOG_LOG, s, 1.0);
    }

    /** Default: LOG_LOG using Constants, with given s and D. */
    public static DsResult estimateDs(double[] eigenvalues, double s, double D) {
        return estimateDs(eigenvalues, FitSpace.LOG_LOG, s, D);
    }

    /** Default: uses Constants for all knobs (s=1, D=1). */
    public static DsResult estimateDs(double[] eigenvalues, FitSpace fitSpace) {
        return estimateDs(eigenvalues, fitSpace, 1.0, 1.0);
    }

    /** Default: uses Constants for all knobs, with s (D=1). */
    public static DsResult estimateDs(double[] eigenvalues, FitSpace fitSpace, double s) {
        return estimateDs(eigenvalues, fitSpace, s, 1.0);
    }

    /** Default: uses Constants for all knobs, with s and D. */
    public static DsResult estimateDs(double[] eigenvalues, FitSpace fitSpace, double s, double D) {
        Objects.requireNonNull(eigenvalues, "eigenvalues null");
        Objects.requireNonNull(fitSpace, "fitSpace null");

        return estimateDs(
                eigenvalues,
                fitSpace,
                s,
                D,
                Constants.EIGEN_ZERO_CLAMP,
                Constants.DS_MERGE_DUPLICATES,
                Constants.DS_MERGE_ABS_TOL,
                Constants.DS_MERGE_REL_TOL,
                Constants.DS_MIN_WINDOW_POINTS,
                Constants.DS_USE_RELATIVE_SPAN,
                Constants.DS_MIN_XSPAN_FRACTION,
                Constants.DS_MIN_XSPAN_ABS,
                Constants.DS_TRIM_LO_FRACTION,
                Constants.DS_TRIM_HI_FRACTION
        );
    }

    /**
     * Expert overload + (s, D).
     *
     * For LOG_LOG:
     *   x = log(D) + s*log(lambda)
     * For other fit spaces:
     *   we leave x as lambda (i.e., s and D are ignored) to avoid changing semantics silently.
     *   If you WANT s,D to apply to LINEAR_* fits too, say so and I’ll add it explicitly.
     */
    public static DsResult estimateDs(
            double[] eigenvalues,
            FitSpace fitSpace,
            double s,
            double D,
            double eigenvalueZeroClamp,
            boolean mergeNearDuplicateEigenvalues,
            double mergeAbsTolerance,
            double mergeRelTolerance,
            int minWindowPoints,
            boolean useRelativeSpan,
            double minXSpanFraction,
            double minXSpanAbsolute,
            double trimLoFraction,
            double trimHiFraction
    ) {
        Objects.requireNonNull(eigenvalues, "eigenvalues null");
        Objects.requireNonNull(fitSpace, "fitSpace null");

        if (!Double.isFinite(s) || s <= 0.0) {
            throw new IllegalArgumentException("s must be finite and > 0 (got " + s + ")");
        }
        if (!Double.isFinite(D) || D <= 0.0) {
            throw new IllegalArgumentException("D must be finite and > 0 (got " + D + ")");
        }

        if (eigenvalues.length < 4) {
            throw new IllegalArgumentException("Too few eigenvalues: " + eigenvalues.length);
        }
        requireFiniteNonNegative("eigenvalueZeroClamp", eigenvalueZeroClamp);
        requireFiniteNonNegative("mergeAbsTolerance", mergeAbsTolerance);
        requireFiniteNonNegative("mergeRelTolerance", mergeRelTolerance);
        if (minWindowPoints < 4) {
            throw new IllegalArgumentException("minWindowPoints must be >= 4 (got " + minWindowPoints + ")");
        }
        requireFiniteNonNegative("minXSpanFraction", minXSpanFraction);
        requireFiniteNonNegative("minXSpanAbsolute", minXSpanAbsolute);

        trimLoFraction = clamp01(trimLoFraction);
        trimHiFraction = clamp01(trimHiFraction);
        if (trimLoFraction + trimHiFraction >= 0.9) {
            throw new IllegalArgumentException("trim fractions too large: lo+hi must be < 0.9");
        }

        final int nTotal = eigenvalues.length;

        // 1) Filter finite entries, split into "zero-like" and "positive-for-fit"
        double[] positives = new double[nTotal];
        int posCount = 0;
        int zeroLikeCount = 0;

        for (double v : eigenvalues) {
            if (!Double.isFinite(v)) continue;

            if (Math.abs(v) <= eigenvalueZeroClamp) {
                zeroLikeCount++;
            } else if (v > eigenvalueZeroClamp) {
                positives[posCount++] = v;
            } else {
                // negative: ignore (Laplacian numerical noise)
            }
        }

        if (posCount < 4) {
            throw new IllegalArgumentException(
                    "Not enough positive eigenvalues > clamp (" + eigenvalueZeroClamp + "). count=" + posCount
            );
        }

        double[] lambda = Arrays.copyOf(positives, posCount);
        Arrays.sort(lambda);

        // 2) Unique + multiplicities (optional merge)
        final double maxLam = lambda[lambda.length - 1];
        final double tol = mergeNearDuplicateEigenvalues
                ? Math.max(mergeAbsTolerance, mergeRelTolerance * Math.max(1.0, maxLam))
                : 0.0;

        UniqueSpectrum spec = (mergeNearDuplicateEigenvalues)
                ? uniqueWithCountsMerged(lambda, tol)
                : uniqueWithCountsExact(lambda);

        if (spec.unique.length < 4) {
            throw new IllegalArgumentException("Spectrum too small after unique/merge: uniqueCount=" + spec.unique.length);
        }

        // 3) Build (x,y) points from UNIQUE eigenvalues; rho from cumulative multiplicity + zeroLikeCount.
        int U = spec.unique.length;
        double[] x = new double[U];
        double[] y = new double[U];

        final double logD = (fitSpace == FitSpace.LOG_LOG) ? Math.log(D) : 0.0;

        int cumulative = zeroLikeCount;
        for (int j = 0; j < U; j++) {
            double lam = spec.unique[j];
            int mult = spec.counts[j];
            cumulative += mult;

            double rho = (double) cumulative / (double) nTotal;

            switch (fitSpace) {
                case LOG_LOG -> {
                    double logLam = Math.log(lam);
                    x[j] = logD + s * logLam;
                    y[j] = Math.log(rho);
                }
                case LINEAR_LOG -> {
                    x[j] = lam;              // s,D intentionally not applied here
                    y[j] = Math.log(rho);
                }
                case LINEAR_LINEAR -> {
                    x[j] = lam;              // s,D intentionally not applied here
                    y[j] = rho;
                }
                default -> throw new IllegalStateException("Unknown fitSpace: " + fitSpace);
            }
        }

        // 4) Trimming
        int L = U;
        int loTrim = (int) Math.floor(trimLoFraction * L);
        int hiTrim = (int) Math.floor(trimHiFraction * L);
        int scanStart = Math.min(Math.max(0, loTrim), L);
        int scanEndExclusive = Math.max(scanStart, L - hiTrim);

        int Lscan = scanEndExclusive - scanStart;
        if (Lscan < 4) {
            throw new IllegalArgumentException("Too few points after trimming: Lscan=" + Lscan +
                    " (L=" + L + ", trimLo=" + trimLoFraction + ", trimHi=" + trimHiFraction + ")");
        }

        int minPts = Math.min(minWindowPoints, Lscan);
        if (minPts < 4) {
            throw new IllegalArgumentException("minWindowPoints too large for available points: requested=" +
                    minWindowPoints + " available=" + Lscan);
        }

        // 5) Choose minXSpan
        // IMPORTANT: for LOG_LOG, x spans scale by factor s, so we scale the absolute minimum span too.
        double minXSpanAbsEff = (fitSpace == FitSpace.LOG_LOG) ? (s * minXSpanAbsolute) : minXSpanAbsolute;

        double minXSpan = chooseMinXSpan(
                x,
                scanStart,
                scanEndExclusive,
                fitSpace,
                useRelativeSpan,
                minXSpanFraction,
                minXSpanAbsEff
        );

        // 6) Best-window scan
        Best best = scanBestWindowByAdjR2WithPenalty(x, y, scanStart, scanEndExclusive, minPts, minXSpan);

        if (best == null) {
            throw new IllegalStateException(
                    "Failed to estimate: no window satisfied minPts=" + minPts + " and minXSpan=" + minXSpan +
                            " (uniquePoints=" + Lscan + ", fitSpace=" + fitSpace +
                            ", trim=[" + trimLoFraction + "," + trimHiFraction + "], s=" + s + ", D=" + D + ")"
            );
        }

        double ds = (fitSpace == FitSpace.LOG_LOG) ? (2.0 * best.fit.slope) : Double.NaN;

        return new DsResult(
                ds,
                best.fit.adjR2,
                best.fit.slope,
                best.fit.intercept,
                fitSpace,
                best.start,
                best.endExclusive
        );
    }

    public static long[] quantize(double[] eigenvalues) {
        Objects.requireNonNull(eigenvalues, "eigenvalues null");

        final double scale = Constants.EIGEN_DEFAULT_SCALE;
        final double clamp = Constants.EIGEN_ZERO_CLAMP;

        long[] q = new long[eigenvalues.length];
        for (int i = 0; i < eigenvalues.length; i++) {
            double v = eigenvalues[i];
            if (!Double.isFinite(v)) {
                q[i] = 0L;
                continue;
            }
            if (Math.abs(v) < clamp) v = 0.0;
            double scaled = v * scale;
            if (scaled > Long.MAX_VALUE || scaled < Long.MIN_VALUE) {
                throw new ArithmeticException("Quantization overflow at i=" + i + " value=" + v);
            }
            q[i] = Math.round(scaled);
        }
        return q;
    }

    public static double[] dequantize(long[] quantizedEigenvalues) {
        Objects.requireNonNull(quantizedEigenvalues, "quantizedEigenvalues null");
        final double scale = Constants.EIGEN_DEFAULT_SCALE;

        double[] out = new double[quantizedEigenvalues.length];
        for (int i = 0; i < out.length; i++) out[i] = quantizedEigenvalues[i] / scale;
        return out;
    }

    private static final class UniqueSpectrum {
        final double[] unique;
        final int[] counts;

        UniqueSpectrum(double[] unique, int[] counts) {
            this.unique = unique;
            this.counts = counts;
        }
    }

    private static UniqueSpectrum uniqueWithCountsExact(double[] sorted) {
        if (sorted.length == 0) return new UniqueSpectrum(new double[0], new int[0]);

        double[] u = new double[sorted.length];
        int[] c = new int[sorted.length];
        int m = 0;

        double cur = sorted[0];
        int cnt = 1;

        for (int i = 1; i < sorted.length; i++) {
            double x = sorted[i];
            if (x == cur) {
                cnt++;
            } else {
                u[m] = cur;
                c[m] = cnt;
                m++;
                cur = x;
                cnt = 1;
            }
        }
        u[m] = cur;
        c[m] = cnt;
        m++;

        return new UniqueSpectrum(Arrays.copyOf(u, m), Arrays.copyOf(c, m));
    }

    private static UniqueSpectrum uniqueWithCountsMerged(double[] sorted, double tol) {
        if (sorted.length == 0) return new UniqueSpectrum(new double[0], new int[0]);
        requireFiniteNonNegative("tol", tol);

        double[] u = new double[sorted.length];
        int[] c = new int[sorted.length];
        int m = 0;

        int i = 0;
        while (i < sorted.length) {
            double sum = sorted[i];
            int cnt = 1;
            int j = i + 1;

            while (j < sorted.length && Math.abs(sorted[j] - sorted[j - 1]) <= tol) {
                sum += sorted[j];
                cnt++;
                j++;
            }

            u[m] = sum / (double) cnt;
            c[m] = cnt;
            m++;

            i = j;
        }

        return new UniqueSpectrum(Arrays.copyOf(u, m), Arrays.copyOf(c, m));
    }

    public static final class Best {
        final int start;
        final int endExclusive;
        final Fit fit;
        final double score;

        Best(int start, int endExclusive, Fit fit, double score) {
            this.start = start;
            this.endExclusive = endExclusive;
            this.fit = fit;
            this.score = score;
        }
    }

    public static final class Fit {
        final double slope;
        final double intercept;
        final double r2;
        final double adjR2;

        Fit(double slope, double intercept, double r2, double adjR2) {
            this.slope = slope;
            this.intercept = intercept;
            this.r2 = r2;
            this.adjR2 = adjR2;
        }
    }

    private static Best scanBestWindowByAdjR2WithPenalty(double[] x, double[] y, int lo, int hi, int minPts, double minXSpan) {
        final double tieEps = 1e-12;
        final double penalty = 0.01;

        Best best = null;
        double bestScore = Double.NEGATIVE_INFINITY;
        double bestAdj = Double.NEGATIVE_INFINITY;
        int bestLen = -1;
        int bestStart = Integer.MAX_VALUE;

        for (int s = lo; s <= hi - minPts; s++) {
            for (int e = s + minPts; e <= hi; e++) {

                double span = x[e - 1] - x[s];
                if (!(span >= minXSpan)) continue;

                Fit fit = fitLineKahan(x, y, s, e);
                if (!Double.isFinite(fit.adjR2)) continue;

                int len = e - s;
                double score = fit.adjR2 - (penalty / (double) len);

                boolean better = false;

                if (score > bestScore + tieEps) {
                    better = true;
                } else if (Math.abs(score - bestScore) <= tieEps) {
                    if (fit.adjR2 > bestAdj + tieEps) {
                        better = true;
                    } else if (Math.abs(fit.adjR2 - bestAdj) <= tieEps) {
                        if (len > bestLen) {
                            better = true;
                        } else if (len == bestLen) {
                            if (s < bestStart) better = true;
                        }
                    }
                }

                if (better) {
                    bestScore = score;
                    bestAdj = fit.adjR2;
                    bestLen = len;
                    bestStart = s;
                    best = new Best(s, e, fit, score);
                }
            }
        }

        return best;
    }

    private static Fit fitLineKahan(double[] x, double[] y, int a, int b) {
        int n = b - a;
        if (n < 2) return new Fit(Double.NaN, Double.NaN, Double.NaN, Double.NaN);

        double sx = 0.0, cx = 0.0;
        double sy = 0.0, cy = 0.0;

        for (int i = a; i < b; i++) {
            double vx = x[i];
            double yx = vx - cx;
            double tx = sx + yx;
            cx = (tx - sx) - yx;
            sx = tx;

            double vy = y[i];
            double yy = vy - cy;
            double ty = sy + yy;
            cy = (ty - sy) - yy;
            sy = ty;
        }

        double mx = sx / n;
        double my = sy / n;

        double sxx = 0.0, cxx = 0.0;
        double sxy = 0.0, cxy = 0.0;
        double syy = 0.0, cyy = 0.0;

        for (int i = a; i < b; i++) {
            double dx = x[i] - mx;
            double dy = y[i] - my;

            double vxx = dx * dx;
            double vxy = dx * dy;
            double vyy = dy * dy;

            double k1 = vxx - cxx;
            double t1 = sxx + k1;
            cxx = (t1 - sxx) - k1;
            sxx = t1;

            double k2 = vxy - cxy;
            double t2 = sxy + k2;
            cxy = (t2 - sxy) - k2;
            sxy = t2;

            double k3 = vyy - cyy;
            double t3 = syy + k3;
            cyy = (t3 - syy) - k3;
            syy = t3;
        }

        if (!(sxx > 0.0) || !(syy > 0.0)) {
            return new Fit(Double.NaN, Double.NaN, Double.NaN, Double.NaN);
        }

        double slope = sxy / sxx;
        double intercept = my - slope * mx;

        double ssRes = 0.0, cRes = 0.0;
        for (int i = a; i < b; i++) {
            double yhat = slope * x[i] + intercept;
            double r = y[i] - yhat;
            double vr = r * r;

            double k = vr - cRes;
            double t = ssRes + k;
            cRes = (t - ssRes) - k;
            ssRes = t;
        }

        double r2 = 1.0 - (ssRes / syy);

        if (Double.isFinite(r2)) {
            if (r2 < -1e-12) r2 = -1e-12;
            if (r2 > 1.0 + 1e-12) r2 = 1.0;
        }

        int p = 1;
        double adjR2;
        if (n <= p + 1) {
            adjR2 = Double.NaN;
        } else {
            adjR2 = 1.0 - (1.0 - r2) * (n - 1.0) / (n - p - 1.0);
        }

        return new Fit(slope, intercept, r2, adjR2);
    }

    private static double chooseMinXSpan(
            double[] x,
            int scanStart,
            int scanEndExclusive,
            FitSpace fitSpace,
            boolean useRelativeSpan,
            double minXSpanFraction,
            double minXSpanAbsolute
    ) {
        if (scanEndExclusive - scanStart < 2) return minXSpanAbsolute;

        double xmin = x[scanStart];
        double xmax = x[scanStart];
        for (int i = scanStart + 1; i < scanEndExclusive; i++) {
            double xi = x[i];
            if (xi < xmin) xmin = xi;
            if (xi > xmax) xmax = xi;
        }
        double totalSpan = xmax - xmin;

        if (fitSpace == FitSpace.LOG_LOG && useRelativeSpan && Double.isFinite(totalSpan) && totalSpan > 0.0) {
            double rel = minXSpanFraction * totalSpan;
            return Math.max(minXSpanAbsolute, rel);
        }
        return minXSpanAbsolute;
    }

    private static void requireFiniteNonNegative(String name, double v) {
        if (!Double.isFinite(v) || v < 0.0) {
            throw new IllegalArgumentException(name + " must be finite and >= 0 (got " + v + ")");
        }
    }

    private static double clamp01(double x) {
        if (!Double.isFinite(x)) return 0.0;
        if (x < 0.0) return 0.0;
        if (x > 1.0) return 1.0;
        return x;
    }
}
