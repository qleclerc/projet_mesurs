# projet_mesurs
Code repository for the mesurs group project

## model folder

Contains the script for the model.

### Equations of the model

We modelled SARS-CoV-2 transmission in a company of $N$ employees using a compartmental model. In this model, employees can be susceptible to the respiratory disease $S$, exposed to the disease but not yet infectious $E$, infectious and asymptomatic $I_A$, infectious and symptomatic $I_S$, or recovered $R$. 


$$\frac{dS}{dt} = - \lambda S I_A$$

$$\frac{dE}{dt} = \lambda S I_A - \sigma E$$

$$\frac{dI_A}{dt} = p_A \sigma E - \gamma_A I_A$$

$$\frac{dI_S}{dt} = (1-p_A) \sigma E - \gamma_A I_S$$

$$\frac{dR}{dt} = \gamma_A I_A + \gamma_S I_S$$


In this compartmental model, infected individuals can develop symptoms with probability $1-p_A$ but their incubation period (here, time from infection to time of the onset of infectiousness) is the same whether they develop symptoms or not and is equal to $\frac{1}{\sigma}$. Infectious and symptomatic individuals are infectious for $\frac{1}{\gamma_S}$ days and are assumed to be on medical leave from their symptom/infectiousness onset to their recovery. Thereby, these individuals do not contribute to the propagation of the epidemic within the company. Infectious and asymptomatic individuals are infectious for $\frac{1}{\gamma_A}$ days and are responsible for the spread of the disease within the company.

The transmission rate $\lambda$ is divided into three terms as follows:

$$\lambda = \frac{5}{7} (1-\alpha) \beta \frac{I_A}{N} + \frac{5}{7} \alpha \epsilon \lambda_v + \frac{2}{7} \lambda_v$$

Where $\beta$ is expressed using SARS-CoV-2 $R_0$ that we derived from the next-generation matrix for a frequency-dependent model $\beta = \frac{R_0 \gamma_A}{p_A}$, $\alpha$ is the proportion of employees teleworking, $\lambda_v$ is the transmission rate from the community, and $\epsilon$ is a coefficient reducing the transmission from the community on teleworking days.

In addition to the transmission process of the infectious disease, we modelled the number of individuals who will ultimately develop a chronic disease following exposure to teleworking. To do so, we stratified the compartmental model into two populations, one population that will not develop a chronic disease, and a second that will develop a chronic disease. Given the different time scales of occurrence of the chronic disease and the infectious disease, we assume that the two populations mix homogeneously between them.

### Parameterization of the model

According to [Wu et al., 2022](https://doi.org/10.1001/jamanetworkopen.2022.28008), the pooled incubation period for all SARS-CoV-2 variants is 6.57 days. In this case, $\sigma = \frac{1}{6.57}$. 

| SARS-CoV-2 variant | Mean incubation period (in days) |
| :----------------: | :------------------------------: |
| Alpha              | 5.00                             |
| Beta               | 4.50                             |
| Delta              | 4.41                             |
| Omicron            | 3.42                             |

[Koelle et al., 2022](https://doi.org/10.1126/science.abm4915) have reviewed how our understanding of the transmission of SARS-CoV-2 has changed over the pandemic mostly using modelling studies. The first studies in Wuhan estimated that the $R_0$ lied between 2 and 4 before the lockdown. In their review, [Dhungel et al., 2022]( https://doi.org/10.3390/ijerph191811613) have estimated a pooled $R_0$ of 2.66 from studies published in the early months of the pandemic. 

[Byambasuren et al., 2020](https://doi.org/10.3138/jammi-2020-0030) estimated in the first year of the pandemic that the proportion of asymptomatic cases does not exceed 20% and is on average equal to 17%. This would lead to $p_A = 0.17$

Need to find estimates for $\gamma_A$ and $\gamma_S$. However, is 

[Li et al., 2020](https://doi.org/10.1126/science.abb3221): asymptomatic individuals are 55% as infectious as symptomatic individuals which should reduce the $R_0$ value.

## analysis folder

Contains scripts used to run the model.

## figures folder

Contains figures produced by analysis scripts.

