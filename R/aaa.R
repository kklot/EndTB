pars_moldova <- c(
    N          = 23842802, # population
    beta_s     = 4.4,      # Mean rate of transmission per TB case with DS-TB: est. 4.4 (95% CrI 2.4-6)
    beta_r     = 12,       # 
    kappa      = 0.1,      # U(0.1 – 10) Relative infectivity, post vs pre- care-seeking
    b          = 0.01,     # Birth rate: WHO (GHO) (15) – adjusted to yield 1% annual population growth from 1970
    mu         = 0.015,    # Background mortality rate: WHO (GHO)(15), corresponds to mean life expectancy of 68 years
    mu_tb      = 0.15,     # 0,15 (0.14-0.18)
    theta_s    = 0.14,     # Vynnycky E et al. (7)
    theta_r    = 0.14,     # Proportion infections being ‘fast’ progressors to active disease
    rho        = 0.001,    # Horsburgh CR et al. (8) Breakdown to active disease
    sigma      = 7.9,      # est - 7.9 (95% CrI 3.7-11.8) Per-capita hazard of developing symptoms
    delta      = 7.5,      # (95% CrI 2.7-11.5), Per-capita rate of initial presentation to care
    gamma      = 12,       # 12 (95% CrI 8.5-15) Interval between care-seeking episodes
    phi        = 52,       # Treatment initiation delay
    varepsilon = 0.93,     # Probability of diagnosis per patient-provider interaction
    omega      = 0.93,     # Treatment initiation probability, WHO country report 2016-assumption to match TB coverage of 87%, assuming that TB coverage is the product of omega and varepsilon
    tau_0      = 2,        # Treatment duration, first line,
    tau_1      = 0.5,      # Treatment duration, second line,
    chi_s      = 0.2,      # Treatment default rate 
    chi_r      = 0.25,     # Moldova’s National Tb programme 2016-2017
    varrho     = 0.05,     # 0.05 (95% CrI 0.005-0.09) MDR acquisition rate while on FL
    # Relapse, per-capita hazard rates
    r_0        = 0.032,    # relapse following treatment completion
    r_1        = 0.14,     # relapse following treatment default
    r_2        = 0.0015,   # relapse >2 years after treatment 
    r_3        = 0.7,      # relapse in MDR recovered from 1st line
    varsigma   = 0.5,      # ‘Stabilisation’ of relapse risk following treatment, Based in Thomas et al (10): most relapse occurs in first 2 yr after treatment.
    # Treatment success rate (Moldova’s National TB programme 2016-2017)
        c_s0   = 0.8,    # DS-TB in first line treatment, 
        c_r0   = 0.3,      # MDR-TB in first line treatment 
        c_r1   = 0.48,     # MDR-TB in second line treatment
    # Drug sensitivity testing rates (DST) - ECDC 2016 (28)
    m_n        = 0.74,     # new cases
    m_r        = 0.47,     # retreated cases
    xi         = 0.8       # Second line treatment initiation probability
)

pars_moldova_null <- c(
    N          = 23842802, # population
    beta_s     = 4.4,      # Mean rate of transmission per TB case with DS-TB: est. 4.4 (95% CrI 2.4-6)
    beta_r     = 0,        # 
    kappa      = 1,        # U(0.1 – 10) Relative infectivity, post vs pre- care-seeking
    b          = 0.01,     # Birth rate: WHO (GHO) (15) – adjusted to yield 1% annual population growth from 1970
    mu         = 0.01,     # Background mortality rate: WHO (GHO)(15), corresponds to mean life expectancy of 68 years
    mu_tb      = 0.15,     # 0,15 (0.14-0.18)
    theta_s    = 0.14,     # Vynnycky E et al. (7)
    theta_r    = 0.14,     # Proportion infections being ‘fast’ progressors to active disease
    rho        = 0.001,    # Horsburgh CR et al. (8) Breakdown to active disease
    sigma      = 7.9,      # est - 7.9 (95% CrI 3.7-11.8) Per-capita hazard of developing symptoms
    delta      = 0,        # (95% CrI 2.7-11.5), Per-capita rate of initial presentation to care
    gamma      = 0,        # 12 (95% CrI 8.5-15) Interval between care-seeking episodes
    phi        = 0,        # Treatment initiation delay
    varepsilon = 0,        # Probability of diagnosis per patient-provider interaction
    omega      = 0,        # Treatment initiation probability, WHO country report 2016-assumption to match TB coverage of 87%, assuming that TB coverage is the product of omega and varepsilon
    tau_0      = 0,        # Treatment duration, first line,
    tau_1      = 0,        # Treatment duration, second line,
    chi_s      = 0,        # Treatment default rate 
    chi_r      = 0,        # Moldova’s National Tb programme 2016-2017
    varrho     = 0,        # 0.05 (95% CrI 0.005-0.09) MDR acquisition rate while on FL
    # Relapse, per-capita hazard rates
    r_0        = 0,        # relapse following treatment completion
    r_1        = 0,        # relapse following treatment default
    r_2        = 0,        # relapse >2 years after treatment 
    r_3        = 0,        # relapse in MDR recovered from 1st line
    varsigma   = 0,        # ‘Stabilisation’ of relapse risk following treatment, Based in Thomas et al (10): most relapse occurs in first 2 yr after treatment.
    # Treatment success rate (Moldova’s National TB programme 2016-2017)
        c_s0   = 0,        # DS-TB in first line treatment, 
        c_r0   = 0,        # MDR-TB in first line treatment 
        c_r1   = 0,        # MDR-TB in second line treatment
    # Drug sensitivity testing rates (DST) - ECDC 2016 (28)
    m_n        = 0,        # new cases
    m_r        = 0,        # retreated cases
    xi         = 0         # Second line treatment initiation probability
)

init_moldova <- c(
    S   = 23842802,
    LS1 = 1,
    AS1 = 0,
    IS1 = 0,
    DS1 = 0,
    ES1 = 0,
    FS1 = 0,
    AS2 = 0,
    IS2 = 0,
    DS2 = 0,
    ES2 = 0,
    FS2 = 0,
    RS0 = 0,
    RS1 = 0,
    RS2 = 0,
    LR1 = 0, # MDR
    AR1 = 0,
    IR1 = 0,
    DR1 = 0,
    ER1 = 0,
    FR1 = 0,
    AR2 = 0,
    IR2 = 0,
    DR2 = 0,
    ER2 = 0,
    FR2 = 0,
    RR0 = 0,
    RR1 = 0,
    RR2 = 0,
    G1 = 0,
    G2 = 0, 
    NOT = 0,
    MOR = 0,
    FOIs = 0,
    FOIr = 0,
    PREV = 0
)

states_moldova <- c(
    S   = "Susceptible",
    LS1 = "Sen-1st-Latent",
    AS1 = "Sen-1st-Asymptomatic",
    IS1 = "Sen-1st-Symptomatic",
    DS1 = "Sen-1st-Presented",
    ES1 = "Sen-1st-Seekcare",
    FS1 = "Sen-1st-Treat",
    AS2 = "Sen-2nd-Asymptomatic",
    IS2 = "Sen-2nd-Symptomatic",
    DS2 = "Sen-2nd-Presented",
    ES2 = "Sen-2nd-Seekcare",
    FS2 = "Sen-2nd-Treat",
    RS0 = "Sen-2nd-Recover-Stable",
    RS1 = "Sen-2nd-Relapse-Low",
    RS2 = "Sen-2nd-Relapse-High",
    LR1 = "Res-1st-Latent",
    AR1 = "Res-1st-Asymptomatic",
    IR1 = "Res-1st-Symptomatic",
    DR1 = "Res-1st-Presented",
    ER1 = "Res-1st-Seekcare",
    FR1 = "Res-1st-Treat",
    AR2 = "Res-2nd-Asymptomatic",
    IR2 = "Res-2nd-Symptomatic",
    DR2 = "Res-2nd-Presented",
    ER2 = "Res-2nd-Seekcare",
    FR2 = "Res-2nd-Treat",
    RR0 = "Res-2nd-Recover-Stable",
    RR1 = "Res-2nd-Relapse-Low",
    RR2 = "Res-2nd-Relapse-High",
    G1  = "Res-1st-1st",
    G2  = "Res-1st-2nd",
    NOT = "Notifications",
    MOR = "Mortality",
    FOIs = "FOI DS",
    FOIr = "FOI DR",
    PREV = "Prevalence"
)