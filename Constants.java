package iso;

public class Constants 
{
    public static double EIGEN_DEFAULT_SCALE = 1e6; //1e10
    public static double EIGEN_ZERO_CLAMP = 0.5 / EIGEN_DEFAULT_SCALE;
    
    public static double SPECTRAL_EXPONENT = 1.0;
    public static double DIFFUSION_CONSTANT = 1.0;
    
    public static double DS_EPS = 1e-6;   //1e10
    public static boolean DS_MERGE_DUPLICATES = true;
    public static double DS_MERGE_ABS_TOL = EIGEN_ZERO_CLAMP;
    public static double DS_MERGE_REL_TOL = 1e-10; // or 1e-12
    public static int DS_MIN_WINDOW_POINTS = 16;
    public static boolean DS_USE_RELATIVE_SPAN = true;
    public static double DS_MIN_XSPAN_FRACTION = 0.20;
    public static double DS_MIN_XSPAN_ABS = 0.2;   // key change (was 1.0)
    public static double DS_TRIM_LO_FRACTION = 0.01;
    public static double DS_TRIM_HI_FRACTION = 0.6;   
    
    public static long CURV_COEFF_QUANT_BUCKET = 128;
    /** Controls whether coefficient fitting is done in log-scaled space (recommended). */
    public static boolean CURV_FIT_IN_LOG_SPACE = true;
    /** Minimum K used before log() to avoid -Inf when CURV_FIT_IN_LOG_SPACE=true. */
    public static double CURV_LOG_FIT_MIN_K = 1e-300;
    /** If true AND SamplingDomain=LOG_T, time grid uses Chebyshev nodes in log(t) for best conditioning. */
    public static boolean CURV_TGRID_USE_CHEBYSHEV_NODES_LOGT = true;
    /** Clip for exp/log to keep numbers sane. */
    public static double CURV_LOG_EXP_CLIP = 700.0;
    /** Clip threshold for exp(-t*mu): if -t*mu < -CURV_EXP_CLIP treat as zero contribution. */
    public static double CURV_EXP_CLIP = 700.0;
    /** Safe bounds for exp(logMu) when computing effective mu = D*lambda^s via logs. */
    public static double CURV_LOG_MU_MIN = -745.0;
    public static double CURV_LOG_MU_MAX = 709.0;
    /** Degeneracy grouping tolerances in effective spectrum mu-space. */
    public static double CURV_MU_GROUP_ABS_TOL = EIGEN_ZERO_CLAMP;
    public static double CURV_MU_GROUP_REL_TOL = 1e-10; //1e-12
    public static double CURV_MU_GROUP_ULP_MULT = 16.0;
    /** Time grid construction tolerances/limits (effective spectrum mu-space). */
    public static double CURV_TGRID_EIG_TOL = 1e-9; //1e-9
    public static double CURV_TGRID_SPAN_MAX = 1e6;
    public static double CURV_TGRID_MAX_LAM_T = 1000.0;
    public static double CURV_TGRID_MIN_RATIO = 1.05;
    /** Auto time-grid sizing: T = pointsPerCoeff * p, clamped. */
    public static int CURV_TGRID_POINTS_PER_COEFF = 6;
    public static int CURV_TGRID_MIN_POINTS = 32;
    public static int CURV_TGRID_MAX_POINTS = 512;
    
    /** Default quantization scale for curvature/heat coefficients. */
    public static double CURV_COEFF_DEFAULT_SCALE = 1e6;  //1e10
    
    /** Signature / BFS layering */
    public static int SIGNATURE_BFS_LAYERS = 1000000;
}
