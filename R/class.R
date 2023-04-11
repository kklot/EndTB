#' Tuberculosis Models
#' 
#' @description TB model
#' @export
TBM <- R6::R6Class("TBM",
list(
    #' @field model TBM in TMB code
    model = NULL,
    #' @field current_pars the parameters used for model simulation
    current_pars = pars_moldova_null,
    #' @description Init the model
    #' @param pars A named vectors of parameters, must matched with the existing parameters (see `pars_moldova` for details)
    #' @param tmax Maximum number of years to simulate
    #' @param dt Time steps, default to 1
    #' @return NULL
    initialize = function(pars = NULL, tmax = 150, dt = 1) {
        if (!is.null(pars)) self$update_pars(pars)
        self$model <- TMB::MakeADFun(
            data = list(init = init_moldova, tmax = tmax, dt = dt),
            parameters = list(pars = self$current_pars),
            DLL = "EndTB", silent = T
        )
    },
    #' @description update model's parameters
    #' @param pars a named vector of parameters 
    update_pars = function(pars) {
        self$current_pars <<- replace(self$current_pars, names(pars), pars)
    },
    #' @description plot model dynamic
    #' @param states list of model states to plot
    #' @return a `ggplot` object.
    plot = function(states = c("S", "AS1", "LS1", "IS1")) {
        requireNamespace("ggplot2", quietly = TRUE)
        p <- private$out_long() %>%
        dplyr::filter(name %in% states) %>%
        ggplot2::ggplot(ggplot2::aes(time, value)) +
        ggplot2::geom_line(ggplot2::aes(color = state)) +
        ggplot2::facet_wrap(~state, nrow = 2, scales = "free") +
        ggplot2::labs(x = "Year", y = "%") +
        ggplot2::theme(legend.position = 'none')
        print(p)
        invisible(p)
    }
),
list(
    out_long = function() {
        self$output %>%
            t() %>%
            dplyr::as_tibble(.name_repair = function(x) make.names(x, unique = TRUE)) %>%
            dplyr::rename_with(~ c("time", names(init_moldova))) %>%
            dplyr::slice(-dplyr::n()) %>%
            dplyr::mutate(dplyr::across(-time, function(x) x / rowSums(dplyr::across(-time)))) %>%
            tidyr::pivot_longer(-time) %>%
            dplyr::mutate(state = factor(name, levels = names(states_moldova), labels = states_moldova))
    }
),
list(
    #' @field output return model's simulation
    output = function() {self$model$report()$out}
)
)