---
title: Notes
markup: "mmark"
output: html_document
---

## {.tabset .tabset-fade .tabset-pills}


### Demonstration  {.tabset .tabset-fade .tabset-pills}

https://www.symbolab.com/

The ASA's Statement on p-Values: Context, Process, and Purpose

http://dx.doi.org/10.1080/00031305.2016.1154108

stplanr

https://www.rdocumentation.org/packages/stplanr




Normal_PDF|Normal_CDF|Chi^2_PDF|Chi^2_CDF
---|---|---|---
![Normal_Distribution_PDF](https://upload.wikimedia.org/wikipedia/commons/thumb/1/1b/Normal_distribution_pdf.png/320px-Normal_distribution_pdf.png) | ![Normal_Distribution_CDF](https://upload.wikimedia.org/wikipedia/commons/thumb/1/19/Normal_distribution_cdf.png/320px-Normal_distribution_cdf.png) | ![Chi-square_distributionPDF](https://upload.wikimedia.org/wikipedia/commons/thumb/2/21/Chi-square_distributionPDF.png/320px-Chi-square_distributionPDF.png) | ![Chi-square_distributionCDF](https://upload.wikimedia.org/wikipedia/commons/thumb/c/cb/Chi-square_distributionCDF.png/320px-Chi-square_distributionCDF.png)
F_PDF|F_CDF|Cauchy_PDF|Cauchy_CDF
![F-distribution_pdf](https://upload.wikimedia.org/wikipedia/commons/thumb/7/74/F-distribution_pdf.svg/320px-F-distribution_pdf.svg.png)|![F-distribution_cdf](https://upload.wikimedia.org/wikipedia/commons/thumb/d/df/F_distributionCDF.png/320px-F_distributionCDF.png)|![Cauchy_pdf](https://upload.wikimedia.org/wikipedia/commons/thumb/b/b7/Cauchy_distribution_pdf.png/320px-Cauchy_distribution_pdf.png)|![Cauchy_cdf](https://upload.wikimedia.org/wikipedia/commons/thumb/e/e4/Cauchy_distribution_cdf.png/320px-Cauchy_distribution_cdf.png) 
Exponential_pmf|Exponential_cdf|Binomial_pmf|Binomial_cdf
![Exponential_distribution_pdf](https://upload.wikimedia.org/wikipedia/commons/thumb/b/b1/Exponential_distribution_pdf.png/320px-Exponential_distribution_pdf.png)|![Exponential_distribution_cdf](https://upload.wikimedia.org/wikipedia/commons/thumb/7/77/Exponential_distribution_cdf.png/320px-Exponential_distribution_cdf.png)|![Binomial_distribution_pmf](https://upload.wikimedia.org/wikipedia/commons/thumb/7/75/Binomial_distribution_pmf.svg/320px-Binomial_distribution_pmf.svg.png)|![Binomial_distribution_cdf](https://upload.wikimedia.org/wikipedia/commons/thumb/5/55/Binomial_distribution_cdf.svg/320px-Binomial_distribution_cdf.svg.png)






#### interactive visualization


- [Interpreting Confidence Intervals](http://rpsychologist.com/d3/CI/) 

- [Interpreting Cohen's d effect size](http://rpsychologist.com/d3/cohend/)

- [Interpreting Correlations](http://rpsychologist.com/d3/correlation/)

- [Understanding Statistical Power and Significance Testing](http://rpsychologist.com/d3/NHST/)

### Univariable Distribution

1. [Chart of Univariable Distribution](https://www.wolfram.com/mathematica/new-in-8/parametric-probability-distributions/index.html)


1. [Eighty Univariate Distributions and Their Relationships Displayed in a Matrix Format](https://ieeexplore-ieee-org.proxy.lib.pdx.edu/abstract/document/5755180)

[Probability Theory and Mathematical Statistics](https://onlinecourses.science.psu.edu/stat414/node/109/)

### Glossary of Statistical Terms


- ANOVA: 
Analysis of variance usually refers to an analysis of a continuous dependent variable where all the predictor variables are categorical. One-way ANOVA, where there is only one predictor variable (factor; grouping variable), is a generalization of the 2-sample t-test. ANOVA with 2 groups is identical to the t-test. Two-way ANOVA refers to two predictors, and if the two are allowed to interact in the model, two-way ANOVA involves cross-classification of observations simultaneously by both factors. It is not appropriate to refer to repeated measures within subjects as two-way ANOVA (e.g., treatment× time). An ANOVA table sometimes refers to statistics for more complex models, where explained variation from partial and total effects are displayed and continuous variables may be included.

- Bayes’ rule or theorem: 

$Pr(A|B) = \frac{Pr(B|A) Pr(A)}{Pr(B)}$, read as the probability that event A happens given that event B has happened equals the probability that B happens given that A has happened multiplied by the (unconditional) probability that A happens and divided by the (unconditional) probability that B happens. Bayes’ rule follows immediately from the law of conditional probability which states that $Pr(A|B) = Pr(A \cap B)Pr(B)$.

- Bayesian inference: 

A branch of statistics based on Bayes’ theorem. Bayesian inference doesn’t use P -values and generally does not test hypotheses. It requires one to formally specify a probability distribution encapsulating the prior knowledge about, say, a treatment effect. The state of prior knowledge can be specified as “no knowledge” by using a flat distribution. Once the prior distribution is specified, the data are used to modify the prior state of knowledge to obtain the post-experiment state of knowledge. Final probabilities computed in the Bayesian framework are probabilities of various treatment effects.

- bias: 
A systematic error. Examples: a miscalibrated machine that reports cholesterol too high by 20mg% on the average; a satisfaction questionnaire that leads patients to never report that they are dissatisfied with their medical care; using each patient’s lowest blood pressure over 24 hours to describe a drug’s antihyptertensive properties.

- binary variable: 

A variable whose only two possible values, usually zero and one.

- bootstrap: 

A simulation technique for studying properties of statistics without the need to have the infinite population available. The most common use of the bootstrap involves taking random samples (with replacement) from the original dataset and studying how some quantity of interest varies. Each random sample has the same number of observations as the original dataset. Some of the original subjects may be omitted from the random sample and some may be sampled more than once. The bootstrap can be used to compute standard deviations and confidence limits without assuming a model. For example, if one took 200 samples with replacement from the original dataset, computed the sample median from each sample, and then computed the sample standard deviation of the 200 medians, the result would be a good estimate of the true standard deviation of the original sample median. The bootstrap can also be used to internally validate a predictive model without holding back patient data during model development.

- calibration: 

Reliability of predicted values, i.e., extent to which predicted values agree with observed values. For a predictive model a calibration curve is constructed by relating predicted to observed values in some smooth manner. The calibration curve is judged against a 45◦ line. Miscalibration could be called bias. Calibration error is frequently assessed for predicted event probabilities. If for example 0.4 of the time it rained when the predicted probability of rain was 0.4, the rain forecast is perfectly calibrated.

- case-control study: 

A study in which subjects are selected on the basis of their outcomes, and then exposures (treatments) are ascertained. For example, to assess the association between race and operative mortality one might select all patients who died after open heart surgery in a given year and then select an equal number of patients who survived, matching on several variables other than race so as to equalize (control for) their distributions between the cases and non-cases.

- categorical variable: 

A variable having only certain possible values for which there is no logical ordering of the values. Also called a nominal, polytomous, discrete categorical variable or factor.

- censoring: 

When the response variable is the time until an event, subjects not followed long enough for the event to have occurred have their event times censored at the time of last follow-up. This kind of censoring is right censoring. For example, in a follow-up study, patients entering the study during its last year will be followed a maximum of 1 year, so they will have their time until event censored at 1 year or less. Left censoring means that the time to the event is known to be less than some value. In interval censoring the time is known to be in a specified interval. Most statistical analyses assume that what causes a subject to be censored is independent of what would cause her to have an event. If this is not the case, informative censoring is said to be present. For example, if a subject is pulled off of a drug because of a treatment failure, the censoring time is indirectly reflecting a bad clinical outcome and the resulting analysis will be biased.

- conditional probability: 

The probability of the veracity of a statement or of an event A given that a specific condition B holds or that an event B has already occurred, denoted by P (A|B). This is a probability in the presence of knowledge captured by B. For example, if the condition B is that a person is male, the conditional probability is the probability of A for males. It could be argued that there is no such thing as a completely unconditional probability. In this example one is implicitly conditioning on humans even if not considering the person’s sex.

- confidence limits: 

To say that the 0.95 confidence limits for an unknown quantity are [a, b] means that 0.95 of similarly constructed confidence limits in repeated samples from the same population would contain the unknown quantity. Very loosely speaking one could say that she is 0.95 “confident” that the unknown value is in the interval [a, b], although in the frequentist school unknown parameters are constants, so they are either inside or outside intervals and there are no probabilities associated with these events. The interpretation of a single confidence interval in frequentist statistics is highly problematic. Note that a confidence interval should be symmetric about a point estimate only when the distribution of the point estimate is symmetric. Many confidence intervals are asymmetric, e.g., intervals for probabilities, odds ratios, and other ratios.

- confounder: 

A variable which is correlated with the response variable and with the treatment assignment (or exposure variable). A confounder, when properly controlled for, can explain away an apparent association between the treatment and the response.

- continuous variable: 

A variable that can take on any number of possible values. Practically speaking, when a variable can take on at least, say, 10 values, it can be treated as a continuous variable. For example, it can be plotted on a scatterplot and certain meaningful calculations can be made using the variable.


- critical value: 

The value of a test statistic (e.g., t, F , χ2, z) that if exceeded by the observed test statistic would result in statistical significance at a chosen α level or better. For a z-test (normal deviate test) the critical level of z is 1.96 when α = 0.05 for a two-sided test. For t and F tests, critical values decrease as the sample size increases, as one requires less penalty for having to estimate the population variance as n gets large.

- cross-validation: 

This technique involves leaving out m patients at a time, fitting a model on the remaining n − m patients, and obtaining an unbiased evaluation of predictive accuracy on the m patients. The estimates are averaged over ≥ n/m repetitions. Cross-validation provides estimates that have more variation than those from bootstrapping. It may require > 200 model fits to yield precise estimates of predictive accuracy.

- detectable difference: 

The value of a true population treatment effect (difference between two treatments) that if held would result in a statistical test having exactly the desired power.

- discrimination: 

A variable or model’s discrimination ability is its ability to separate subjects having a low responses from subjects having high responses. One way to quantify discrimination is the ROC curve area.

- dummy variable: 

A device used in a multivariable regression model to describe a categorical predictor without assuming a numeric scoring. Indicator variable might be a better term. For example, treat- ments A, B, C might be described by the two dummy predictor variables X1 and X2, where X1 is a binary variable taking on the value of 1 if the treatment for the subject is B and 0 otherwise, and X2 takes on the value 1 if the subject is under treatment C and 0 otherwise. The two dummy variables completely define 3 categories, because when X1 = X2 = 0 the treatment is A.

- estimate: 

A statistical estimate of a parameter based on the data. See parameter. Examples include the sample mean, sample median, and estimated regression coefficients.

- frequentist statistical inference: 

Currently the most commonly used statistical philosophy. It uses hy- pothesis testing, type I and II errors, power, P -values, confidence limits, and adjustments of P -values for testing multiple hypotheses from the same study. Final probabilities computed using frequentist methods, P -values, are probabilities of obtaining values of statistics. The frequentist approach is also called the sampling approach as it considers the distribution of statistics over hypothetical repeated samples from the same population. The frequentist approach is concerned with long-run operating char- acteristics of statistics and estimates. Because of this and because of the backwards time/information ordering of P -values, frequentist testing requires complex multiplicity adjustments but provides no guiding principles for exactly how those adjustments should be derived. Frequentist statistics involves massive confusion of two ideas: (1) the apriori probability that an experiment will generate misleading information (e.g., the chance of a false positive result or type I error) and (2) the evidence for an asser- tion after the experiment is run. The latter should not involve a multiplicity adjustment, but because the former does, frequentists do not know how to interpret the latter when multiple hypotheses are tested or when a single hypothesis is tested sequentially.

- goodness of fit: 

Assessment of the agreement of the data with either a hypothesized pattern (e.g., independence of row and column factors in a contingency table or the form of a regression relationship) or a hypothesized distribution (e.g., comparing a histogram with expected frequencies from the normal distribution).

- Hawthorne effect: 

A change in a subject response that results from the subject knowing she is being observed.


- inter-quartile range: 

The range between the outer quartiles (25th and 75th percentiles).

- least squares estimate: 

The value of a regression coefficient that results in the minimum sum of squared errors, where an error is defined as the difference between an observed and a predicted dependent variable value.

- linear regression model: 

This is also called OLS or ordinary least squares and refers to regression for a continuous dependent variable, and usually to the case where the residuals are assumed to be Gaussian. The linear model is sometimes called general linear model, not to be confused with generalized linear model where the distribution can take on many non-Gaussian forms.

- logistic regression model: 

A multivariable regression model relating one or more predictor variables to the probabilities of various outcomes. The most commonly used logistic model is the binary logistic model7,6 which predicts the probability of an event as a function of several variables. There are several types of ordinal logistic models for predicting an ordinal outcome variable, and there is a polytomous logistic model for categorical responses. The binary and polytomous models generalize the χ2 test for testing for association between categorical variables. One commonly used ordinal model, the proportional odds model1, generalizes the Wilcoxon 2-sample rank test. Binary logistic models are useful for predicting events in which time is not very important. They can be used to predict events by a specified time, but this can result in a loss of information. Logistic models are used to estimate adjusted odds ratios as well as probabilities of events.

- maximum likelihood estimate: 

An estimate of a statistical parameter (such as a regression coefficient, mean, variance, or standard deviation) that is the value of that parameter making the data most likely to have been observed. MLEs have excellent statistical properties in general, such as converging to population values as the sample size increases, and having the best precision from among all such competing estimators. When the data are normally distributed, maximum likelihood estimates of regression coefficients and means are equivalent to least squares estimates. When the data are not normally distributed (e.g. binary outcomes, or survival times), maximum likelihood is the standard method to estimate the regression coefficients (e.g. logistic regression, Cox regression).

- mean: 

Arithmetic average, i.e., the sum of all the values divided by the number of observations. The mean of a binary variable is equal to the proportion of ones because the sum of all the zero and one values equals the number of ones. The mean can be heavily influenced by outliers.

- median: 

Value such that half of the observations’ values are less than and half are greater than that value. The median is also called the 50th percentile or the 0.5 quantile. The median is not heavily influenced by outliers so it can be more representative of “typical” subjects. When the data happen to be normally (Gaussian) distributed, the median is not as precise as the mean in describing the central tendency.

- multiple comparisons: 

It is common for one study to involve the calculation of more than one P -value.
For example, the investigator may wish to test for treatment effects in 3 groups defined by disease etiology, she may test the effects on 4 different patient response variables, or she may look for a significant difference in blood pressure at each of 24 hourly measurements. When multiple statistical tests are done, the chances of at least one of them resulting in a false positive finding increases as the number of tests increase. This is called “inflation of type I error.” When one wishes to control the overall type I error, individual tests can be done using a more stringent α level, or individual P -values can be adjusted upward. Such adjustments are usually dictated when using frequentist statistics, as P -values mean the probability of getting a result this impressive if there is really no effect, and “this impressive” can be taken to mean “this impressive given the large number of statistics examined.”

- multivariable model: 
A model relating multiple predictor variables (risk factors, treatments, etc.) to a single response or dependent variable. The predictor variables may be continuous, binary, or cat- egorical. When a continuous variable is used, a linearity assumption is made unless the variable is expanded to include nonlinear terms. Categorical variables are modeled using dummy variables so as to not assume numeric assignments to categories.

- multivariate model: 

A model that simultaneously predicts more than one dependent variable, e.g. a model to predict systolic and diastolic blood pressure or a model to predict systolic blood pressure 5 min. and 60 min. after drug administration.

- nominal significance level: 

In the context of multiple comparisons involving multiple statistical tests, the apparent significance level α of each test is called the nominal significance level. The overall type I error for the study, the probability of at least one false positive result, will be greater than α.

- nonparametric estimator: 

A method for estimating a parameter without assuming an underlying dis- tribution for the data. Examples include sample quantiles, the empirical cumulative distribution, and the Kaplan-Meier survival curve estimator.

- nonparametric tests: 

A test that makes minimal assumptions about the distribution of the data or about certain parameters of a statistical model. Nonparametric tests for ordinal or continuous variables are typically based on the ranks of the data values. Such tests are unaffected by any one-one transformation of the data, e.g., by taking logs. Even if the data come from a normal distribution, rank tests lose very little efficiency (they have a relative efficiency of 3/π = 0.955 if the distribution is normal) compared with parametric tests such as the t-test and the linear correlation test. If the data are not normal, a rank test can be much more efficient than the corresponding parametric test. For these reasons, it is not very fruitful to test data for normality and then to decide between the parametric and nonparametric approaches. In addition, tests of normality are not always very powerful. Examples of nonparametric tests are the 2-sample Wilcoxon-Mann-Whitney test, the 1-sample Wilcoxon signed-rank test (usually used for paired data), and the Spearman, Kendall, or Somers rank correlation tests. Even though nonparametric tests do not assume a specific distribution for a group, they assume a connection between the distributions of any two groups. For example, the logrank test assumes proportional hazards, i.e., that the survival curve for group A is a power of the survival curve for group B. The Wilcoxon test, for optimal power, assumes that the cumulative distributions are in proportional odds.

- normal distribution: 

A symmetric, bell-shaped distribution that is most useful for approximating the distribution of statistical estimators. Also called the Gaussian distribution. The normal distribution cannot be relied upon to approximate the distribution of raw data. The normal distribution’s bell shape follows a rigid mathematical equation of the form e−x^2 . For a normal distribution the probability that a measurement will fall within ±1.96 standard deviations of the mean is 0.95.

- null hypothesis: 

Customarily but not necessarily a hypothesis of no effect, e.g., no reduction in mean blood pressure or no correlation between age and blood pressure. The null hypothesis, labeled H0, is often used in the frequentist branch of statistical inference as a “straw person”; classical statistics often assumes what one hopes doesn’t happen (no effect of a treatment) and attempts to gather evidence against that assumption (i.e., tries to reject H0). H0 usually specifies a single point such as 0mmHg reduction in blood pressure, but it can specify an interval, e.g., H0: blood pressure reduction is between -1 and +1 mmHg. “Null hypotheses” can also be e.g. H0: correlation between X and Y is 0.5.

- observational study: 

Study in which no experimental condition (e.g., treatment) is manipulated by the investigator, i.e., randomization is not used.

- odds: 

The probability an event occurs divided by the probability that it doesn’t occur. An event that
occurs 0.90 of the time has 9:1 odds of occurring since 0.9(1−0.9)= 9.

- odds ratio: 

The odds ratio for comparing two groups (A, B) on their probabilities of an outcome occurring is the odds of the event occurring for group A divided by the odds that it occurs for group B. If P(A) and P(B) represent the probability of the outcome for the two groups of subjects, the A : B odds ratio is (PA/1−PA)/(PB/1−PB). Odds ratios are in the interval [0, ∞). An odds ratio for a treatment is a measure of relative effect of that treatment on a binary outcome. As summary measures, odds ratios have advantages over risk ratios: they don’t depend on which of two possible outcomes is labeled the “event”, and any odds ratio can apply to any probability of outcome in the reference group. Because of this, one frequently finds that odds ratios for comparing treatments are relatively constant across different types of patients. The same is not true of risk ratios or risk differences; these depend on the level of risk in the reference group.

- one-sided test: 

A test designed to test a directional hypothesis, yielding a one-sided P -value. For example, one might test the null hypothesis H0 that there is no difference in mortality between two treatments, with the alternative hypothesis being that the new drug lowers mortality. See also two-sided test.

- ordinal variable: 

A categorical variable for which there is a definite ordering of the categories. For example, severity of lower back pain could be ordered as none, mild, moderate, severe, and coded using these names or using numeric codes such as 0,1,2,10. Spacings between codes are not important.

- P -value: 

The probability of getting a result (e.g., t or χ2 statistics) as or more extreme than the observed statistic had H0 been true. An α-level test would reject H0 if P ≤ α. However, the P -value can be reported instead of choosing an arbitrary value of α. Examples:   
(1) An investigator compared two randomized groups for differences in systolic blood pressure, with the two mean pressures being 134.4 mmHg and 138.2 mmHg. She obtained a two-tailed P = 0.03. This means that if there is truly no difference in the population means, one would expect to find a difference in means exceeding 3.8 mmHg in absolute value 0.03 of the time. The investigator might conclude there is evidence for a treatment effect on mean systolic blood pressure if the statistical test’s assumptions are true.  
(2) An investigator obtained P = 0.23 for testing a correlation being zero, with the sample correlation being 0.08. The probability of getting a correlation this large or larger in absolute value if the population correlation is zero is 0.08. No conclusion is possible other than (a) more data are needed and (b) there is no convincing evidence for or against a zero correlation. For both of these examples confidence intervals would be helpful.

- paired data: 

When each subject has two response measurements, there is a natural pairing to the data and the two responses are correlated. The correlation results from the fact that generally there is more variation between subjects than there is within subjects. Sometimes one can take the difference or log ratio of the two responses for each subject, and then analyze these “effect measures” using an unpaired one-sample approach such as the Wilcoxon signed-rank test or the paired t-test. One must be careful that the effect measure is properly chosen so that it is independent of the baseline value.

- parameter: 

An unknown quantity such as the population mean, population variance, difference in two means, or regression coefficient.

- parametric model: 

A model based on a mathematical function having a few unknown parameters.

- parametric test: 

A test which makes specific assumptions about the distribution of the data or specific assumptions about model parameters. Examples include the t-test and the Pearson product-moment linear correlation test.

- percentile: 

The p-th percentile is the value such that np/100 of the observations’ values are less than that value. The p-th quantile is the value such that np of the observations’ values are less.

- posterior probability: 

In a Bayesian context, this is the probability of an event after making use of the information in the data. In other words, it is the prior probability of an event after updating it with the data. Posterior probability can also be called post-test probability if one equates a diagnostic test with “data” (see also ROC curve).

- power: 

Probability of rejecting the null hypothesis for a set value of the unknown effect. Power could also be called the sensitivity of the statistical test in detecting that effect. Power increases when the sample size and true unknown effect increase and when the inter-subject variability decreases. In a two-group comparison, power generally increases as the allocation ratio gets closer to 1:1. For a given experiment it is desirable to use a statistical test expected to have maximum power (sensitivity). A less powerful statistical test will have the same power as a better test that was applied after discarding some of the observations. For example, testing for differences in the proportion of patients with hypertension in a 500-patient study may yield the same power as a 350-patient study which used blood pressure as a continuous variable.

- precision: 

Degree of absence of random error. The precision of a statistical estimator is related to the expected error that occurs when approximating the infinite-data value. In other words, when you try to estimate some measure in a population, the precision is related to the error in the estimate. So precision can be thought of as a “margin of error” in estimating some unknown value. Precision can be quantified by the width of a confidence interval and sometimes by a standard deviation of the estimator (standard error). For the confidence intervals, a “margin for error” is computed so that the quoted interval has a certain probability of containing the true value (e.g., population mean difference). Some authors define precision as the reciprocal of the variance of an estimate. By that definition, precision increases linearly as the sample size increases. If instead one defines precision on the original scale of measurement instead of its square (i.e., if one uses the standard error or width of a confidence interval), precision increases as the square root of the sample size.

- predictor, explanatory variable, risk factor, covariate, covariable, independent variable: 

quan- tities which may be associated with better or worse outcome.

- prior probability: 

The probability of an event as it could best be assessed before the experiment. In diagnostic testing this is called the pre-test probability. The prior probability can come from an objective model based on previously available information, or it can be based on expert opinion. In some Bayesian analyses, prior probabilities are expressed as probability distributions which are flat lines, to reflect a complete absence of knowledge about an event. Such distributions are called non- informative, flat, or reference distributions, and analyses based on them fully let the data “speak for themselves.”

- probability: 

In the frequentist school, the probability of an event denotes the limit of the long-term fraction of occurrences of the event. This notion of probability implies that the same experiment which generated the outcome of interest can be repeated infinitely often. Even a coin will change after 100,000 flips. Likewise, some may argue that a patient is “one of a kind” and that repetitions of the same experiment are not possible. One could reasonably argue that a “repetition” does not denote the same patient at the same stage of the disease, but rather any patient with the same severity of disease (measured with current technology). There are other schools of probability that do not require the notion of replication at all. For example, the school of subjective probability (associated with the Bayesian school) “considers probability as a measure of the degree of belief of a given subject in the occurrence of an event or, more generally, in the veracity of a given assertion” (see P. 55 of5). de Finetti defined subjective probability in terms of wagers and odds in betting. A risk-neutral individual would be willing to wager $P that an event will occur when the payoff is $1 and her subjective probability is P for the event. The domain of application of probability is all-important. We assume that the true event status (e.g., dead/alive) is unknown, and we also assume that the information the probability is conditional upon (e.g. Pr{death | male, age=70}) is what we would check the probability against. In other words, we do not ask whether Pr(death | male, age=70) is accurate when compared against Pr(death | male, age=70, meanbp=45, patient on downhill course).

- prospective study: 

One in which the study is first designed, then the subjects are enrolled. Prospective studies are usually characterized by intentional data collection.

- quartiles: 

The 25th and 75th percentiles and the median. The three values divide a variables distributions into four intervals containing equal numbers of observations.

- random error: 

An error caused by sampling from a group rather than knowing the true value of a quantity such as the mean blood pressure for the entire group, e.g., healthy men over age 80. One can also speak of random errors in single measurements for individual subjects, e.g., the error in using a single blood pressure measurement to represent a subject’s long-term blood pressure.

- random sample: 

A sample selected by a random device that ensures that the sample (if large enough) is representative of the infinite group. A probability sample is a kind of random sample in which each possible subject has a known probability of being sampled, but the probabilities can vary. For example, one may wish to over-sample African-Americans in a study to ensure good representation. In that case one could sample African-Americans with probability of 1.0 and others with a probability of 0.5.

- randomness: 

Absence of a systematic pattern. One might wish to examine whether some hormone level varies systematically over the day as opposed to having a random pattern, or whether events such as epileptic seizures tend to cluster or occur randomly in time. Sometimes the residuals in an ordinary regression model are plotted against the order in which subjects were accrued to make sure that the pattern is random (e.g., there was no learning trend for the investigators).

- rate: 

A ratio such as a change per unit time. Rates are often limits, and shouldn’t be confused with probabilities. The latter are constrained to be between 0 and 1 whereas there are no constraints on possible values for rates.

- regression to the mean: 

Tendency for a variable that has an extreme value on its first measurement to have a more typical value on its second measurement. For example, suppose that subjects must have LDL cholesterol > 190mg% to qualify for a study, and the median LDL cholesterol for qualifying subjects at the screening visit was 230 mg%. The median LDL cholesterol value at their second visit might be 200mg%, with several of the subjects having values below 190. This is the “sophomore slump” in baseball; second-year players are watched when they have phenomenal rookie years. Regression to the mean also takes many other forms, all arising because variables or subgroups are not examined at random but rather because they appear “impressive”: (1) One might compare 5 treatments with a control and choose the treatment having the maximum difference. On a repeated study that treatment’s average response will be found to be much closer to that of the control. (2) In a randomized controlled trial the investigators may wish to estimate the effect of treatment in multiple subgroups. They find that in 40 left-handed diabetics the treatment multiplies mortality by 0.4. If the study is replicated, they would find that the mortality reduction in left-handed diabetics is much closer to the mortality reduction in the overall sample of patients. (3) Researchers study the association between 40 possible risk factors and some outcome, and find that the factor with the strongest association had a correlation of 0.5 with the response. On replication, the correlation will be much lower. This result is very related to what happens in stepwise variable selection, where the most statistically significant variables selected will have their importance (regression coefficients) greatly overstated.

- residual: 

A statistical quantity that should be unrelated to certain other variables because their effects should have already been subtracted out. In ordinary multiple regression, the most commonly used residual is the difference between predicted and observed values.

- retrospective study: 

A study in which subjects were already enrolled before the study was designed, or the outcome of interest has occurred before the start of the study (an in a case control study). Such studies often have difficulties such as absence of needed adjustment (confounder) variables and missing data.

- risk: 

Often used as another name for probability but a more accurate definition is the probability of an adverse event × the severity of the loss that experiencing that event would entail.

- semi-parametric: 

‘Parametric’ assumptions may be made about some aspects of a model, while other components may be estimated ‘non-parametrically’. In the Cox regression procedure, a parametric model for the relative hazard is overlaid on a non-parametric estimate of baseline hazard2.

- significance level: 

A preset value of α against which P -values are judged in order to reject H0 (see Type I error). Sometimes a P -value itself is called the significance level.

- standard deviation: 

A measure of the variability (spread) of measurements across subjects. The standard deviation has a simple interpretation only if the data distribution is Gaussian (normal), and in that restrictive case the mean ±1.96 standard deviations is expected to cover 0.95 of the distribution of the measurement. Standard deviation is the square root of the variance.

- standard error: 

The standard deviation of a statistical estimator. For example, the standard deviation of a mean is called the standard error of the mean, and it equals the standard deviation of individual measurements divided by the square root of the sample size. Standard errors describe the precision of a statistical summary, not the variability across subjects. Standard errors go to zero as the sample size → ∞.

- survival analysis: 

A branch of statistics dealing with the analysis of the time until an event such as death. Survival analysis is distinguished by its emphasis on estimating the time course of events and in dealing with censoring. See Cox model.

- survival function: 

The probability of being free of the event at a specified time.

- survival time: 

Interval between the time origin and the occurrence of the event or censoring2.

- symmetric distribution: 

One in which values to the left of the mean by a certain amount are just as likely to be observed as values to the right of the mean by the same amount. For symmetric distributions, the population mean and median are identical and the distance between the 25th and 50th percentiles equals the distance between the 50th and 75th percentiles.

- two-sided test: 

A test that is non-directional and that leads to a two-sided P -value. If the null hypothesis H0 is that two treatments have the same mortality outcome, a two-sided alternative is that the mortality difference is nonzero. Two-sided P -values are larger than one-sided P -values (they are double if the distribution of the test statistic is symmetric). They can be thought of as a multiplicity adjustment that would allow a claim to be made that a treatment lowers or raises mortality. See also one-sided test.

- type I error: 

False positive rate – the probability of rejecting H0 (i.e., declaring “statistical significance”) when the null hypothesis is in fact true. The type I error is often called α.

- type II error: 

Failing to detect an effect that is real, i.e., the false negative rate. The type II error is referred to as β, which is one minus the power of the test. In other words, the power of the test is 1 − β.

- variance: 

A measure of the spread or variability of a distribution, equaling the average value of the squared difference between measurements and the population mean measurement. From a sample of measurements, the variance is estimated by the sample variance, which is the sum of squared differences from the sample mean, divided by the number of measurements minus 1. The minus 1 is a kind of “penalty” that corrects for estimating the population mean with the sample mean. Variances are typically only useful when the measurements follow a normal or at least a symmetric distribution.



### Machine Learning by PwC

<img src="https://usblogs.pwc.com/emerging-technology/wp-content/uploads/2017/05/machine-learning-overview-thumb.png" style="max-width:15%;min-width:40px;" >
[Machine learning overview](http://usblogs.pwc.com/emerging-technology/machine-learning-overview-infographic-redirect/)

<img src="https://usblogs.pwc.com/emerging-technology/wp-content/uploads/2017/05/machine-learning-methods-thumb.png" style="max-width:15%;min-width:40px;" >
[Machine learning methods](http://usblogs.pwc.com/emerging-technology/machine-learning-methods-infographic/)

<img src="https://usblogs.pwc.com/emerging-technology/wp-content/uploads/2017/05/machine-learning-evolution-thumb.png" style="max-width:15%;min-width:40px;" >
[Machine learning evolution](https://usblogs.pwc.com/emerging-technology/machine-learning-evolution-infographic/)


