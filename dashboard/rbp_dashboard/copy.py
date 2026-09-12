"""Prose the dashboard shows: term definitions, and the caveats that must travel with figures.

Kept here rather than inline so that a claim appears once. The wording follows the manuscript;
where it compresses, it compresses toward the weaker statement.
"""

GLOSSARY = {
    "AUROC": (
        "Area under the ROC curve. The chance that a randomly chosen bound window scores above "
        "a randomly chosen unbound one. 0.5 is coin-flipping, 1.0 is perfect."
    ),
    "negative set": (
        "The windows called 'not bound'. They are chosen, not observed, and this study is about "
        "how much that choice changes the answer."
    ),
    "composition baseline": (
        "A model given only 19 numbers describing nucleotide composition: how much A, C, G and T, "
        "GC fraction, and dinucleotide frequencies. No sequence order, no learning about motifs."
    ),
    "nested contribution": (
        "AUROC of (composition features + the model's score) minus AUROC of composition alone. "
        "What the model adds once plain composition is already accounted for."
    ),
    "apparent AUROC": (
        "The score a paper would print as its headline: composition features plus the model's "
        "score, with no baseline subtracted."
    ),
    "GC-matched": (
        "Unbound windows picked so their overall G+C fraction resembles the bound windows'."
    ),
    "dinucleotide-matched": (
        "A stricter match: pairs of adjacent bases are matched too, not just single-base content."
    ),
    "bias-aware": (
        "Unbound windows drawn to control the technical biases of the assay itself, rather than "
        "matched on composition."
    ),
    "two-stage estimator": (
        "Fit the model, then measure what it adds. The literature surveyed here does this. It "
        "lets information from the evaluation folds leak into the score being evaluated."
    ),
    "cross-fitted estimator": (
        "The model's score for each window is produced by a model that never saw that window's "
        "fold. This closes the leak. It is the primary estimator in this paper."
    ),
    "estimator floor": (
        "What an estimator reports when the true answer is known to be zero. Anything above zero "
        "is the estimator measuring itself."
    ),
    "protein-clustered bootstrap": (
        "Resampling whole proteins rather than single datasets, so two datasets for the same "
        "protein cannot be treated as independent evidence."
    ),
    "span": (
        "Largest divided by smallest across the three negative-set protocols. A span of 1.0 "
        "would mean the protocol made no difference."
    ),
    "held-out benchmark": (
        "135 datasets from a separate published deposit, sharing no dataset with the 94 analysed "
        "here, used to test whether the finding survives outside this pipeline."
    ),
}

# Every figure comparing model classes must carry this. The 4-mer result is cross-fitted and
# primary; the neural results are two-stage and exploratory, and conflating them would overstate
# what the study establishes.
ESTIMATOR_CAVEAT = (
    "The 4-mer result is reported under both estimators, and its cross-fitted value is the "
    "paper's primary estimand. The CNN and SpliceBERT spans are computed with the conventional "
    "two-stage estimator only and are exploratory. They were not cross-fitted, so they carry the "
    "same outer-fold channel the 2-mer test below shows to be non-zero."
)

PANEL_CAVEAT = (
    "95 datasets were selected for the panel. 94 carry all three negative-set protocols and are "
    "the denominator for every three-protocol comparison here."
)

NOT_ONLY_NEGATIVES = (
    "The negative-set construction is what varies, but it is not literally the only thing that "
    "changes. The model and the composition baseline are refitted inside each protocol, and each "
    "protocol keeps the positives its matcher could pair, so retained positive subsets differ "
    "slightly between arms. The manuscript quantifies that difference and shows it does not "
    "account for the result."
)

REPRO_CAVEAT = (
    "This is verification, not complete raw-to-result reproduction. The release checks each "
    "committed table against its recorded evidence and provenance, and every table carries a "
    "checksum. How far a given table can actually be regenerated varies by table: some rebuild "
    "from raw inputs, some only from committed intermediates, some came off cloud jobs that "
    "would have to be rerun, and two have no recorded producing script at all. No claim of "
    "bit-for-bit end-to-end reproducibility is made."
)

NOT_MONITORING = (
    "A static reader over committed CSV files. It is not production monitoring, it is not "
    "connected to any cloud service, and it does not run analyses or recompute estimates."
)


# Every view opens with one of these. The dashboard was legible to someone who already knew the
# result and opaque to everyone else, which is the usual failure of a dashboard built by the
# person who ran the analysis. `checks` is what to actually look at, in order.
#
# `filters` records whether the sidebar filter changes anything on the view. Four of seven views
# show only panel-level published estimates, which correctly do not move when a protein is
# deselected. Without this stated, filtering appears broken.
GUIDE = {
    "Overview": {
        "lead": "A model's measured contribution depends on a choice most papers do not "
                "report: which windows were called unbound. This quantifies that dependence "
                "across 94 ENCODE eCLIP datasets.",
        "checks": [
            "The primary span, 4.84x. That is how far the measured contribution moves across the "
            "three negative-set constructions, with the model class, source peaks, fold "
            "design and estimator held fixed.",
            "The chart below: the bias-aware protocol is the LOWEST bar. It is also the one "
            "with the highest apparent AUROC, which is the counter-intuitive part.",
            "Known-null bias removed, 96%. The conventional estimator reports a positive "
            "result on a question whose answer is zero.",
        ],
        "filters": "static",
    },
    "Why this exists": {
        "lead": "Two things sent me here. A survey of what published methods actually report, "
                "and an earlier study of mine where the controls kept beating the model.",
        "checks": [
            "The survey: seven methods, 2014 to 2025. Look at the right-hand column. Not one "
            "reports what plain letter-counting scores on its own data.",
            "Without that floor, a published AUROC of 0.85 cannot be distinguished from a "
            "baseline of 0.84 plus a rounding error.",
            "The ladder below: conservation alone, with no model at all, scores 0.908 against "
            "the model's 0.829.",
            "And a DIFFERENT protein's model still reaches 0.680, where chance is 0.500. Most "
            "of what looked protein-specific was not.",
        ],
        "filters": "static",
    },
    "Does the fix work?": {
        "lead": "The paper recommends rescaling scores by the headroom left above the "
                "composition baseline. I pre-specified how to test whether that works, then "
                "tested it. This is the result, including the part that failed.",
        "checks": [
            "All three point estimates are positive. Rank agreement improves on every protocol "
            "pair, which is why the recommendation is made at all.",
            "Now the intervals. Only ONE of the three clears zero before correction.",
            "After Bonferroni correction for testing three pairs, none of them clears zero.",
            "That is the honest reading: the direction is consistent, the effect is not "
            "established. The paper says so rather than quietly reporting the uncorrected "
            "interval.",
        ],
        "filters": "static",
    },
    "How negatives are built": {
        "lead": "The unbound windows are not observed, they are constructed. This view is the "
                "construction itself: how a bound window gets paired with an unbound one, how "
                "closely each protocol matches, and what that matching costs.",
        "checks": [
            "The diagram: a bound window is measured, a pool is searched, a partner is "
            "kept. The two composition-matched protocols differ in how strict the match "
            "must be.",
            "The bias-aware row is not a third strictness setting. It draws from a "
            "different pool, other proteins' binding sites, and matches fold rather than "
            "composition.",
            "The match-quality curve: GC-matched and dinucleotide-matched rise steeply. "
            "Bias-aware does not, because it is not matching composition at all.",
            "20% versus 95%. That single design difference is why the bias-aware arm has the "
            "highest composition baseline and so the smallest room left for a model to add.",
            "The cost chart: stricter matching makes the model's raw score WORSE on almost "
            "every dataset, while its measured contribution goes up.",
            "The redraw chart: five random draws land inside one standard error, so none of "
            "this rests on one lucky sample.",
        ],
        "filters": "static",
    },
    "Protocol sensitivity": {
        "lead": "Three ways of choosing the unbound windows, applied to the same 94 datasets "
                "with the model class, peaks, folds and estimator held fixed.",
        "checks": [
            "On Levels: the dinucleotide-matched bars are lowest for apparent AUROC. It looks "
            "like the hardest protocol.",
            "Then on the second chart: dinucleotide-matched is HIGHEST for contribution. "
            "Difficulty and contribution move in opposite directions.",
            "On Per dataset: the distributions overlap heavily. The panel means are a real "
            "difference between protocols, not a clean separation of every dataset.",
            "On The inverse relation: three clouds, each a protocol, stepping down to the "
            "right. That is the same finding as a scatter.",
        ],
        "filters": "partial",
    },
    "Model comparison": {
        "lead": "Whether the protocol dependence is a property of one weak model or of the "
                "measurement itself.",
        "checks": [
            "Press PLAY. Watch the bars change height while their ORDER does not. All three "
            "model classes rank the protocols the same way.",
            "The spans: 3.72x, 5.42x, 7.42x. None is close to 1.0, which would mean the "
            "protocol made no difference.",
            "SpliceBERT contributes most in every protocol and still spans 3.72x. A better "
            "model does not remove the problem.",
        ],
        "filters": "partial",
    },
    "Dataset explorer": {
        "lead": "All 94 datasets, filterable and rankable, with a drill-down into any one.",
        "checks": [
            "On Map: every line is one experiment crossing the three protocols. They shift "
            "together, which is the protocol setting the level.",
            "Then read the counter: 59 of 94 follow the panel ordering exactly and 35 do not. "
            "The effect is systematic in the panel mean, not uniform across datasets.",
            "Use the sidebar filters here. This is the view they were built for.",
            "On Ranking: compare the live descriptive mean of your filtered subset against the "
            "published panel estimate printed beneath it.",
            "On Single dataset: pick a small dataset, then a large one. The small ones are "
            "far noisier, which is why the study's unit of inference is the panel.",
        ],
        "filters": "responds",
    },
    "Cross-fitting": {
        "lead": "The estimator tested against a question whose answer is known to be zero. "
                "This is the methodological core of the paper.",
        "checks": [
            "The first chart: both bars should be at zero. The amber one is not. That is the "
            "conventional estimator reporting signal that cannot exist.",
            "The three percentages: 95.7%, 100.7%, 97.6% of that false signal removed by "
            "cross-fitting.",
            "The second chart: on the REAL 4-mer signal, cross-fitting moves the estimate "
            "slightly UP, not down. So it is not a blanket deflation.",
        ],
        "filters": "static",
    },
    "External validation": {
        "lead": "135 datasets from a separate published deposit, sharing no dataset with the "
                "94 analysed here. Pre-fixed pass criteria, tested once.",
        "checks": [
            "The green panel: protocol dependence replicated at 1.69x, criteria met.",
            "The red panel: the inverse relation did NOT replicate. This is the honest half.",
            "The two charts side by side: negative-1 is higher on BOTH. Same direction, where "
            "the internal panel showed opposite directions.",
        ],
        "filters": "static",
    },
    "Reproducibility": {
        "lead": "How each committed table was produced, and the difference between verifying "
                "an artefact and reproducing it from raw data.",
        "checks": [
            "The provenance mix: 25 of 166 tables are raw-reproducible. The rest are "
            "recomputable from committed evidence or frozen.",
            "The two unattributed tables. They are counted rather than hidden.",
            "The caveat at the top. No claim of bit-for-bit end-to-end reproduction is made.",
        ],
        "filters": "static",
    },
}

# Only two states say anything, and both say the same thing. A view where the sidebar filter
# does nothing gets no chip and no sentence: seven of the ten were carrying a notice about a
# thing that was not happening, which is noise on seven views to buy clarity on three.
FILTER_STATE = {
    "responds": ("filters apply here",
                 "Every chart on this view narrows with the sidebar filters."),
    "partial": ("filters apply here",
                "Per-dataset charts narrow with the sidebar filters and are marked LIVE. "
                "Published panel estimates are fixed values over all 94 datasets and do not "
                "move when you filter."),
    "static": (None, None),
}


# One sentence per chart, for a reader who does not work in this field. Written to answer "what
# am I looking at" before "what does it mean", and to survive being read on its own. A test
# asserts every chart the app draws has an entry here.
PLAIN = {
    "protocol_contribution":
        "Each bar is <b>how much the model added</b> beyond simple letter-counting, under "
        "one way of choosing the comparison sequences. Taller means the model helped more. "
        "All three bars use the same model class on the same source peaks, so they ought to "
        "be close. They are not.",
    "protocol_levels":
        "The darker bar is the score you get from <b>counting letters alone</b>. The lighter bar "
        "is the score with the model added. The gap between them is what the model is worth.",
    "contribution_distribution":
        "Every dot is one experiment. The box shows where the middle half of them sit, so you "
        "can see that the three groups <b>overlap heavily</b> even though their averages differ.",
    "paired_slopes":
        "Each thin line is a single experiment measured <b>both ways</b>. Lines sloping up mean "
        "that experiment scored higher under the protocol on the right.",
    "apparent_versus_contribution":
        "Left-to-right is how <b>impressive the headline score</b> looks. Bottom-to-top is how "
        "much the model actually added. If those agreed, the dots would run diagonally upward. "
        "They run the other way.",
    "animated_model_walk":
        "Press play. The bars change height as the <b>type of model</b> changes, but their "
        "left-to-right order stays the same, which means the pattern is not about one model "
        "being weak.",
    "model_spans":
        "How far each model's answer moved when the negative-set construction changed. "
        "<b>1.0 would mean no effect at all.</b> None of them is close to 1.0.",
    "model_by_protocol":
        "The same three models, grouped by protocol. Note that <b>all three are lowest</b> in "
        "the right-hand group, so they agree about which protocol looks least favourable.",
    "dataset_heatmap":
        "One row per experiment, three columns for the three protocols. Lighter means the model "
        "added more. Hard to read at 94 rows, which is why the trajectory view above is the "
        "default; kept here because it shows the whole panel at once.",
    "protocol_trajectories":
        "Every experiment is one line crossing the three protocols. The thick line is the panel "
        "median. Watch the lines <b>shift together</b> as they cross: that common movement is "
        "the protocol setting the level, and it is what the 4.84x span measures.",
    "ranked_datasets":
        "The experiments where the model added the most, under the protocol you selected. "
        "Hover any bar for the protein and cell line.",
    "dataset_profile":
        "A single experiment, broken out by <b>all three models and all three protocols</b>. "
        "One dataset on its own is noisy, so read this as an illustration rather than evidence.",
    "known_null":
        "A deliberate test with a <b>known answer of zero</b>. Both bars should sit on the line. "
        "The one that does not is the conventional method reporting something that cannot exist.",
    "estimator_comparison":
        "The same two methods on a <b>real</b> signal rather than a fake one. Here they nearly "
        "agree, which is why the problem above is about the method and not about the data.",
    "external_arms":
        "The same measurement repeated on <b>135 completely different experiments</b> from "
        "another research group. Each dot is one of them.",
    "external_levels":
        "Read these two charts together. If the pattern found here also held there, the taller "
        "bar on the left would pair with the <b>shorter</b> bar on the right. It does not.",
    "filtered_versus_published":
        "The solid bar is the <b>published result over all 94 experiments</b>, with its "
        "uncertainty range. The striped bar is just the average of whatever you filtered to, "
        "with no uncertainty range, because a filtered average is not a result.",
    "animated_estimator_walk":
        "Press play. The bars shift as the measuring method is corrected, and the ratio in the "
        "corner drops from 5.42 to 4.84. What matters is that the <b>order does not change</b>: "
        "correcting the method changes the size of the effect, not its direction.",
    "animated_protocol_distribution":
        "Press play to step through the three protocols. Each shape is the spread of 94 "
        "experiments. Watch the whole distribution <b>slide sideways</b> as the only thing that "
        "changed is which comparison sequences were used.",
    "match_quality_curve":
        "How close a match each protocol actually achieved. Further left and higher means a "
        "tighter match. The bottom line is the bias-aware protocol, which <b>does not try</b> "
        "to match composition, and that single design choice drives most of what follows.",
    "match_gap_bars":
        "The typical gap in letter composition between a bound window and the unbound one "
        "paired with it. <b>Shorter bars mean a stricter pairing.</b> The faded bar behind each "
        "is the worst tenth of pairs.",
    "animated_matching":
        "The same curve, one protocol at a time. Press play and watch the shape <b>collapse</b> "
        "when it reaches the bias-aware protocol, which is not matching composition at all.",
    "matching_cost":
        "One bar per experiment. Below the line means the model scored <b>worse</b> when the "
        "comparison windows were matched more strictly. Almost all of them are below the line, "
        "even though the measured contribution goes up.",
    "redraw_stability":
        "The comparison windows are picked at random, so the obvious worry is that one lucky "
        "draw produced the result. Five different random draws are shown. They land <b>well "
        "inside</b> the shaded uncertainty band, so the answer is not an artefact of one sample.",
    "survey_timeline":
        "Seven published methods for predicting where proteins bind RNA, from 2014 to 2025. "
        "Each one builds its comparison set somehow. <b>None of them reports what simple "
        "letter-counting scores on the same data</b>, so none of their headline numbers can be "
        "placed against a floor.",
    "variant_ladder":
        "The study this one grew out of. Read top to bottom: a control using <b>no model at "
        "all</b> scores highest, and scoring each mutation with a completely different "
        "protein's model still lands well above chance. The model was adding something real, "
        "and far less than it looked.",
    "recommendation_intervals":
        "Whether the paper's own suggested fix actually works. Each bar is an uncertainty "
        "range. <b>A range touching the vertical line means the improvement could be nothing.</b> "
        "The blue ranges are before correcting for testing three things at once; the gold ones "
        "are after.",
    "class_ratio_robustness":
        "How many unbound examples you pair with each bound one is a free choice, and one-to-one "
        "is just a convention. Here it is varied four ways. <b>The effect does not go away</b>, "
        "and the ranking of the three protocols holds at every balance.",
    "provenance_mix":
        "How each of the 166 saved data files came to exist. More bars toward the top means more "
        "files that can be <b>rebuilt from scratch</b> rather than simply trusted.",
}


# A narrator, in the author's voice, one passage per view. The guide above says what to check;
# this says why it matters and what it cost to find out. Written to be read by someone who does
# not work on RNA, and kept to claims the tables support.
NARRATOR = {
    "Overview": (
        "Proteins stick to RNA at particular places, and a model can learn to predict where. "
        "To say whether the model is any good you need something to compare against: stretches "
        "of RNA where the protein does **not** stick. Nobody observes those. You choose them.\n\n"
        "I set out to test whether one of these models could flag disease-causing mutations. "
        "Then I built the control, and the control kept beating the model. What I could not "
        "settle was how much of that was the model and how much was my own choice of "
        "comparison sequences, because no paper I could find reported what a plain "
        "letter-counter scores on its own data.\n\n"
        "So I stopped chasing a better number and measured the thing the number rests on. "
        "That is what this is."
    ),
    "Why this exists": (
        "I did not set out to write this paper.\n\n"
        "I wanted to know whether a model that learns where an RNA-binding protein sits could "
        "also flag disease-causing mutations. For a while it looked like it could: SpliceBERT "
        "separated pathogenic from benign non-coding variants at 0.83, having never seen a "
        "disease label.\n\n"
        "Then I built the controls. Evolutionary conservation alone, with no model at all, "
        "reached 0.91. And scoring each variant with a **different** protein's model still "
        "reached 0.68, where chance is 0.50. The model was adding something real, but far less "
        "than it looked, and most of what looked protein-specific was not.\n\n"
        "What kept getting in the way while I checked that was the negatives. Every number I "
        "computed rested on windows I had chosen to call 'not bound', and when I went looking "
        "for what a plain nucleotide counter scores on the same data, I could not find a single "
        "paper that reported it. Seven methods, eleven years, none of them.\n\n"
        "So I stopped chasing a better number and measured the thing the number rests on."
    ),
    "Does the fix work?": (
        "It would be a tidier paper if I left this view out.\n\n"
        "Having shown that the protocol moves the answer, the obvious next question is whether "
        "you can correct for it. I proposed rescaling by the headroom left above the "
        "composition baseline, wrote down in advance how I would test whether that worked, and "
        "then ran the test.\n\n"
        "The direction is right on all three protocol pairs. The magnitude is not established: "
        "one of three intervals clears zero before correcting for multiple comparisons, and "
        "none clears it after.\n\n"
        "I am reporting the corrected interval and calling the recommendation unproven, because "
        "the alternative is to quote the uncorrected one and hope nobody checks. A "
        "recommendation that fails its own pre-specified test is still a finding."
    ),
    "How negatives are built": (
        "This is the part of the study that is a procedure rather than a number, and it is "
        "where the whole result comes from.\n\n"
        "To score a model you need examples of where the protein does not bind. Nobody "
        "measures those. You take a window the protein does bind, measure something about it, "
        "and go looking for a window that resembles it but is not bound. What you choose to "
        "match on is the entire question.\n\n"
        "Match on overall letter balance and 95% of your pairs come within five percentage "
        "points. Match on adjacent pairs of letters as well and you halve the dinucleotide "
        "mismatch. Or ignore composition and correct instead for the quirks of the laboratory "
        "assay, and only 20% of pairs end up composition-matched at all.\n\n"
        "None of those is wrong. All three are in the literature. But they hand the model three "
        "different problems, and the model's measured worth follows the problem rather than the "
        "model."
    ),
    "Protocol sensitivity": (
        "Three reasonable people would build the comparison set three different ways. One "
        "matches the overall balance of letters. One matches adjacent pairs of letters, which "
        "is stricter. One ignores letters and instead corrects for the quirks of the laboratory "
        "method.\n\n"
        "All three appear in the literature. And the measured value of the same model "
        "class, on the same source peaks, moves almost fivefold depending on which one you "
        "picked. They are not three settings of one knob: the first two search free genomic "
        "intervals at different strictness, while the third draws from other proteins' "
        "binding sites entirely.\n\n"
        "The part I did not expect: the protocol whose headline score looks **best** is the one "
        "where the model contributes **least**. A reader comparing headline scores across two "
        "papers is comparing their negative-set choices, not their models."
    ),
    "Model comparison": (
        "The obvious objection is that this is a story about one weak model. So I ran it again "
        "with a convolutional network, and again with SpliceBERT, a large pre-trained "
        "transformer.\n\n"
        "The pattern held in all three. SpliceBERT contributes the most in every protocol, and "
        "its answer still moves 3.7-fold. A better model does not make the problem go away, "
        "because the problem is in the measurement rather than in the model.\n\n"
        "One caution I want to be plain about: only the 4-mer result went through the corrected "
        "estimator. The two neural numbers use the conventional one and are exploratory."
    ),
    "Dataset explorer": (
        "This is the raw material, one row per experiment. I put it here because a panel average "
        "can hide almost anything, and you should be able to check whether the pattern is real "
        "or a handful of outliers dragging a mean around.\n\n"
        "Filter it. Sort it. Pick one experiment and look at it on its own. What you will find "
        "is that the columns differ from each other more than the rows do, which is the finding "
        "restated as a texture rather than a number.\n\n"
        "Anything you compute by filtering is descriptive only. It has no confidence interval "
        "and I do not present it as a result."
    ),
    "Cross-fitting": (
        "This is the part I am most confident about, because it is the part where the answer was "
        "known in advance.\n\n"
        "I gave the method a question whose answer has to be zero. The baseline already contains "
        "every possible two-letter frequency, so a model that knows only two-letter frequencies "
        "cannot add anything. Zero, by construction.\n\n"
        "The conventional method returned a positive number instead, in more than ninety of "
        "ninety-four experiments, in all three protocols. That is not a small bias at the "
        "margin. On the strictest protocol it was ninety percent of the reported effect. "
        "Correcting it is what the rest of this paper is built on."
    ),
    "External validation": (
        "A result that only exists in the pipeline that produced it is not yet a result. So I "
        "wrote down pass and fail criteria, committed them, and then tested them once on 135 "
        "experiments from a different research group, sharing no dataset with mine.\n\n"
        "The main finding replicated. The protocol still moved the answer, by 1.69-fold, inside "
        "the criteria I had fixed beforehand.\n\n"
        "The second finding did not. On their data the difficulty and the contribution move in "
        "the same direction, where in mine they move in opposite directions. I am showing you "
        "that failure at the same size as the success, because a page that only reports the half "
        "that worked is advertising."
    ),
    "Reproducibility": (
        "I want to be precise about what is and is not being claimed here, because "
        "'reproducible' is used loosely.\n\n"
        "Every table in this study can be re-derived and checked against the frozen copy, and "
        "every one carries a checksum. That is verification: it proves the artefact has not "
        "drifted.\n\n"
        "It is **not** a claim that you can rebuild everything from raw sequencing reads on your "
        "laptop. Some of it came off cloud jobs, some is recomputable only from committed "
        "intermediates, and two tables have no recorded producing script at all. Those two are "
        "counted in the chart rather than quietly dropped."
    ),
}
