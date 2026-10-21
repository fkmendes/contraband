package contraband.clock;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.branchratemodel.UCRelaxedClockModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import contraband.math.NodeMath;
import org.apache.commons.math.MathException;

import java.util.Arrays;

/**
 * Trait-driven log-linear clock model for BEAST2.
 *
 * Mathematical Model
 *
 * For each branch (parent → child), the substitution rate r is determined
 * by the branch mid-point body-mass trait value M (already supplied in
 * log-space by NodeMath) via:
 *
 *   log(r) = β · log(M) + log(α_i) + log(μ)
 *
 * which is equivalent in real-space to:
 *
 *   r = M^β · α_i · μ
 *
 * where:
 *
 *   β       – allometric slope (slope)
 *   α_i     – per-branch raw rate (mode-dependent, see below)
 *   μ       – mean clock rate (meanRate)
 *   log(M)  – branch mid-point log-body-mass from NodeMath
 *
 * Raw-Rate Modes
 *
 *   rates mode (rawRate input)
 *     Each branch stores its raw-rate value directly.
 *     When inLogSpace=true, supply log(α_i); otherwise supply α_i > 0.
 *     Reproduces the original TraitLinearClock behaviour.
 *
 *   rateCategories mode (rateCategories + distr inputs)
 *     Each branch stores a discrete category index c_i ∈ {0, …, K-1}.
 *     The raw rate is derived from the inverse CDF of distr:
 *       log(α_i) = log(distr.inverseCDF((c_i + 0.5) / K))
 *     Log-category rates are computed lazily and cached.
 *     inverseCDF output is always real-space; Math.log() is applied once at cache-fill time.
 *
 * inLogSpace Semantics
 *
 *   rates mode
 *     true  → rawRate parameter stores log(α_i); used directly (no log call).
 *     false → rawRate parameter stores α_i > 0; Math.log() applied in getRawRateFromRates().
 *   meanRate follows the same convention.
 *
 *   rateCategories mode
 *     inLogSpace only affects how meanRate is interpreted.
 *     Category log-rates are always derived via Math.log(inverseCDF(...)).
 *
 * Performance
 *
 *   getRawRate() always returns log(α_i).
 *   getRateForBranch() calls Math.exp() exactly once per branch — no extra log/exp calls.
 *
 * Sampled-Ancestor (Zero-Length) Branches
 *
 *   Returns a finite positive rate; branch-length zero means zero likelihood contribution.
 *
 * Prerequisites
 *
 *   Internal-node log-body-mass expectations in NodeMath must be computed by a
 *   pruning-likelihood pass before getRateForBranch(Node) is called.
 */
@Description("Trait-driven log-linear clock model: log(r) = β·log(M) + log(α_i) + log(μ), " +
        "where log(M) is the branch mid-point log-body-mass, α_i is the per-branch raw rate " +
        "(rates or rateCategories mode), μ is the mean rate, and β is the allometric slope.")
public class TraitLinearClockModel extends UCRelaxedClockModel {

    // ── Inputs ────────────────────────────────────────────────────────────────

    final public Input<NodeMath> nodeMathInput = new Input<>(
            "nodeMath",
            "NodeMath object that provides log-space trait expectations for every node.",
            Input.Validate.REQUIRED);

    final public Input<RealParameter> slopeInput = new Input<>(
            "slope",
            "Allometric slope β in log(r) = β·log(M) + log(α_i) + log(μ). " +
                    "Always supplied in real-space regardless of inLogSpace.",
            Input.Validate.REQUIRED);

    final public Input<RealParameter> rawRateInput = new Input<>(
            "rawRate",
            "Per-branch raw rate α_i (rates mode). Dimension = nBranches or 1 (broadcast). " +
                    "When inLogSpace=true supply log(α_i); otherwise supply α_i > 0. " +
                    "Mutually exclusive with rateCategories.",
            Input.Validate.OPTIONAL);

    final public Input<RealParameter> meanRateInput = new Input<>(
            "meanRate",
            "Overall mean clock rate μ. " +
                    "When inLogSpace=true supply log(μ); otherwise supply μ > 0. " +
                    "Fix at 1 (log-space: 0) to avoid unidentifiability with sigmasq.",
            Input.Validate.REQUIRED);

    final public Input<Integer> traitIndexInput = new Input<>(
            "traitIndex",
            "0-based column index of the trait to use when NodeMath holds multiple traits. Default: 0.",
            0);

    final public Input<Boolean> inLogSpaceInput = new Input<>(
            "inLogSpace",
            "If true, rawRate and meanRate are supplied as log-values; slope is always real-space. Default: false.",
            false);

    // ── Mode ──────────────────────────────────────────────────────────────────

    enum Mode { rates, categories }
    Mode mode = Mode.categories;

    // ── State ─────────────────────────────────────────────────────────────────

    RealParameter    rawRates;
    IntegerParameter categories;
    public IntegerParameter    getCategories()   { return categories; }

    ParametricDistribution distribution;
    public ParametricDistribution getDistribution() { return distribution; }

    Tree tree;
    int  branchCount;
    private boolean recompute = true;

    double slope              = 0.0;   // β, always real-space
    double meanRateInLogSpace = 0.0;   // log(μ)
    double meanRate           = 1.0;   // μ real-space, kept in sync with meanRateInLogSpace

    private double storedSlope;
    private double storedMeanRateInLogSpace;
    private double storedMeanRate;

    boolean inLogSpace = false;

    int      LATTICE_SIZE_FOR_DISCRETIZED_RATES = 100;
    double[] rates;        // lazy cache: α_i (real-space) per category; NaN = not yet computed
    double[] storedRates;

    private NodeMath nodeMath;
    private int      traitIndex;

    // ═══ CENTERING CHANGE — revert by removing this field ═══════════════════════
    private double meanTipLogTrait = 0.0;
    // ════════════════════════════════════════════════════════════════════════════

    // ── Lifecycle ─────────────────────────────────────────────────────────────

    @Override
    public void initAndValidate() {
        tree         = treeInput.get();
        branchCount  = tree.getNodeCount() - 1;
        categories   = categoryInput.get();
        rawRates     = rawRateInput.get();
        distribution = rateDistInput.get();

        inLogSpace = inLogSpaceInput.get();

        if (rawRates != null && categories != null)
            throw new IllegalArgumentException(
                    "TraitLinearClockModel: provide rawRate OR rateCategories, not both.");
        if (rawRates == null && categories == null)
            throw new IllegalArgumentException(
                    "TraitLinearClockModel: one of rawRate or rateCategories is required.");

        // mode detection (categories is the default when categoryInput is provided)
        if (categories == null) mode = Mode.rates;

        switch (mode) {
            case categories: {
                LATTICE_SIZE_FOR_DISCRETIZED_RATES = numberOfDiscreteRates.get();
                if (LATTICE_SIZE_FOR_DISCRETIZED_RATES <= 0)
                    LATTICE_SIZE_FOR_DISCRETIZED_RATES = branchCount;
                Log.info.println("  TraitLinearClockModel: using "
                        + LATTICE_SIZE_FOR_DISCRETIZED_RATES
                        + " rate categories to approximate rate distribution across branches.");

                categories.setDimension(branchCount);
                Integer[] initialCategories = new Integer[branchCount];
                for (int i = 0; i < branchCount; i++)
                    initialCategories[i] = Randomizer.nextInt(LATTICE_SIZE_FOR_DISCRETIZED_RATES);
                IntegerParameter other = new IntegerParameter(initialCategories);
                categories.assignFromWithoutID(other);
                categories.setLower(0);
                categories.setUpper(LATTICE_SIZE_FOR_DISCRETIZED_RATES - 1);

                // NaN sentinel: entry not yet computed
                rates       = new double[LATTICE_SIZE_FOR_DISCRETIZED_RATES];
                storedRates = new double[LATTICE_SIZE_FOR_DISCRETIZED_RATES];
                Arrays.fill(rates,       Double.NaN); //TODO: is it needed? need double check
                Arrays.fill(storedRates, Double.NaN);
            }
            break;

            case rates: {
                Log.info.println(
                        "  TraitLinearClockModel: using real rates for rate distribution across branches.");
                int dim = rawRates.getDimension();
                if (dim != 1 && dim != branchCount)
                    throw new IllegalArgumentException(
                            "TraitLinearClockModel: rawRate dimension (" + dim +
                                    ") must be 1 (broadcast) or equal to nBranches (" + branchCount + ").");
                if (!inLogSpace) rawRates.setLower(0.0);
            }
        }

        slope = slopeInput.get().getValue();

        double rawMeanRate = meanRateInput.get().getValue();
        meanRateInLogSpace = inLogSpace ? rawMeanRate : Math.log(rawMeanRate);
        meanRate           = inLogSpace ? Math.exp(rawMeanRate) : rawMeanRate;

        nodeMath   = nodeMathInput.get();
        traitIndex = traitIndexInput.get();

        // ═══ CENTERING CHANGE — revert by removing this block ═══════════════════════
        // Reads directly from the trait Alignment — safe at init, no pruning pass needed.
        // Comma-split matches processContinuousAlignment() in PruneLikelihoodProcess.
        // TODO: check the meanTipLogTrait calculation does that include missing tips or not
        {
            Alignment traitAlignment = nodeMath.traitsValuesInput.get();
            double sum = 0.0;
            int count = 0;
            for (String taxon : traitAlignment.getTaxaNames()) {
                String seq = traitAlignment.getSequenceAsString(taxon);
                if (seq != null && !seq.trim().equals("?") && !seq.trim().equals("-")) {
                    String[] parts = seq.split(",");
                    if (traitIndex < parts.length) {
                        sum += Double.parseDouble(parts[traitIndex].trim());
                        count++;
                    }
                }
            }
            meanTipLogTrait = count > 0 ? sum / count : 0.0;
            Log.info.println("  TraitLinearClockModel: meanTipLogTrait = "
                    + meanTipLogTrait + " (computed from " + count + " tips)");
        }
        // ════════════════════════════════════════════════════════════════════════════
    }

    // ── getRateForBranch ──────────────────────────────────────────────────────

    /**
     * Returns the substitution rate for the branch whose child is node.
     * log(r) = slope·log(M_mid) + log(α_i) + log(μ).
     */
    @Override
    public double getRateForBranch(Node node) {
        if (node.isRoot()) return 1.0;

        if (recompute) {
            synchronized (this) {
                prepare();
                recompute = false;
            }
        }

        // LOG-SPACE version (original) — uncomment to revert:
        // ORIGINAL — revert to this line to undo centering:
//         double logRate = slope * getBranchLogTrait(node) + getRawRate(node) + meanRateInLogSpace;
        // ═══ CENTERING CHANGE — revert by restoring ORIGINAL line above ══════════════
//         double logRate = slope * (getBranchLogTrait(node) - meanTipLogTrait)
//                 + getRawRate(node) + meanRateInLogSpace;
//         return Math.exp(logRate);
        // ════════════════════════════════════════════════════════════════════════════

        // REAL-SPACE version — getRawRate() now returns α_i directly (1 exp() per branch):
        double traitFactor = Math.exp(slope * (getBranchLogTrait(node) - meanTipLogTrait));
        return traitFactor * getRawRate(node) * meanRate;
    }

//     ── Private rate helpers ──────────────────────────────────────────────────

    /** Dispatches to mode-specific helper. Returns α_i in real-space. */
    private double getRawRate(Node node) {
        switch (mode) {
            case categories: return getRawRateFromCategory(node);
            case rates: default: return getRawRateFromRates(node);
        }
    }

    /**
     * Returns α_i in real-space from the rawRates parameter.
     * Supports broadcast (dimension=1) and per-branch (dimension=branchCount).
     * When inLogSpace=true, applies Math.exp() to convert log(α_i) → α_i.
     */
    private double getRawRateFromRates(Node node) {
        int nodeNr = node.getNr();
        if (nodeNr == branchCount)
            nodeNr = node.getTree().getRoot().getNr();
        int idx = rawRates.getDimension() == 1 ? 0 : nodeNr;
        double val = rawRates.getValue(idx);
        // LOG-SPACE version (original) — uncomment to revert:
//         return inLogSpace ? val : Math.log(val);
        return inLogSpace ? Math.exp(val) : val;
    }

    /**
     * Returns α_i in real-space from the inverse CDF of the rate distribution.
     * α_i = distr.inverseCDF((category + 0.5) / K).
     * Computed lazily; NaN entries are filled on first access and reused until
     * the distribution is marked dirty (prepare() resets all entries to NaN).
     */
    private double getRawRateFromCategory(Node node) {
        int nodeNr = node.getNr();
        if (nodeNr == branchCount)
            nodeNr = node.getTree().getRoot().getNr();
        int category = categories.getValue(nodeNr);
        if (Double.isNaN(rates[category])) {
            try {
                // LOG-SPACE version (original) — uncomment to revert:
//                 rates[category] = Math.log(
//                         distribution.inverseCumulativeProbability(
//                                 (category + 0.5) / rates.length));
                rates[category] = distribution.inverseCumulativeProbability(
                        (category + 0.5) / rates.length);
            } catch (MathException e) {
                throw new RuntimeException(
                        "TraitLinearClockModel: failed to compute inverse CDF for category "
                                + category, e);
            }
        }
        return rates[category];
    }

    /** Returns log-body-mass at branch midpoint: (log(M_child) + log(M_parent)) / 2. */
    private double getBranchLogTrait(Node node) {
        int nodeNr   = node.getNr();
        int parentNr = node.getParent().getNr();
        double[] childExpectation  = nodeMath.getExpectationForNode(nodeNr);
        double[] parentExpectation = nodeMath.getInternalNodeExpectation(parentNr);
        return (childExpectation[traitIndex] + parentExpectation[traitIndex]) / 2.0;
    }

    /** Resets the lazy category log-rate cache (NaN = not yet computed). */
    private void prepare() {
        if (rates != null) Arrays.fill(rates, Double.NaN);
    }

    // ── MCMC state management ─────────────────────────────────────────────────

    @Override
    protected boolean requiresRecalculation() {
        recompute = false;

        // distribution changed → clear the log-rate cache
        if (rateDistInput.get().isDirtyCalculation()) {
            recompute = true;
            return true;
        }

        if (slopeInput.get().somethingIsDirty()) {
            slope     = slopeInput.get().getValue();
            recompute = true;
        }

        if (meanRateInput.get().somethingIsDirty()) {
            double val = meanRateInput.get().getValue();
            meanRateInLogSpace = inLogSpace ? val : Math.log(val);
            meanRate           = inLogSpace ? Math.exp(val) : val;
            recompute = true;
        }

        if (mode == Mode.rates && rawRateInput.get().somethingIsDirty())
            recompute = true;

        if (mode == Mode.categories
                && categoryInput.get() != null
                && categoryInput.get().somethingIsDirty())
            recompute = true;

        return recompute;
    }

    @Override
    public void store() {
        if (rates != null)
            System.arraycopy(rates, 0, storedRates, 0, rates.length);
        storedSlope              = slope;
        storedMeanRateInLogSpace = meanRateInLogSpace;
        storedMeanRate           = meanRate;
        super.store();
    }

    @Override
    public void restore() {
        if (rates != null) {
            double[] tmp = rates;
            rates        = storedRates;
            storedRates  = tmp;
        }
        slope              = storedSlope;
        meanRateInLogSpace = storedMeanRateInLogSpace;
        meanRate           = storedMeanRate;
        super.restore();
    }
}
