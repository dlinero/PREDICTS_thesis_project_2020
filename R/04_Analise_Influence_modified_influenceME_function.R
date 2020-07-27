influence_modified <- function (model, group = NULL, select = NULL, obs = FALSE, gf = "single", 
          count = FALSE, delete = TRUE, ...) 
{
  fixef <- NA
  rm(fixef)
  if (is.null(group) & !obs) {
    stop("Please specify either the 'group' parameter, or specify 'obs=TRUE'")
  }
  if (!is.null(group) & obs) {
    stop("Either specify the 'group' parameter, or specify 'obs=TRUE', but not both.")
  }
  ifelse(as.character(model@call)[3] == "data.update", 
         data.adapted <- model.frame(model), data.adapted <- get(as.character(model@call)[3]))
  original.no.estex <- which(substr(names(fixef(model)), 1, 
                                    6) != "estex.")
  n.pred <- length(fixef(model)[original.no.estex])
  if ("(weights)" %in% names(data.adapted)) {
    names(data.adapted)[names(data.adapted) == "(weights)"] <- as.character(model@call$weights)
  }
  if ("(offset)" %in% names(data.adapted)) {
    names(data.adapted)[names(data.adapted) == "(offset)"] <- as.character(model@call$offset)
  }
  if (sum(grepl("offset", names(data.adapted))) > 0) {
    names(data.adapted)[grep("offset", names(data.adapted))] <- gsub("offset\\(|\\)", 
                                                                     "", names(data.adapted)[grep("offset", 
                                                                                                  names(data.adapted))])
  }
  if (!obs) {
    ##### This is the modified part
    # before it was:
    # grouping.names <- grouping.levels(model, group)
    # Now it is:
    
    # Create a character vector with the levels of the group
    V1 <- grouping.levels(model, group)
    
    # Transform the character vector to data frame in order to use subset
    grouping.names <- as.data.frame(V1) %>%
      # Remove the levels that we don't want to exclude and that are specified in a character 
      # vector called exclude
      subset(V1 %nin% exclude) %>%
      #drop levels
      droplevels() %>%
      #transform to a character vector again
      pull() %>% as.character()
    n.groups <- length(grouping.names)
  }
  if (obs) {
    n.obs <- nrow(data.adapted)
  }
  or.fixed <- matrix(ncol = n.pred, nrow = 1, data = fixef(model)[original.no.estex])
  dimnames(or.fixed) <- list(NULL, names(fixef(model))[original.no.estex])
  or.se <- matrix(ncol = n.pred, nrow = 1, data = se.fixef(model)[original.no.estex])
  dimnames(or.se) <- list(NULL, names(fixef(model))[original.no.estex])
  or.vcov <- as.matrix(vcov(model)[original.no.estex, original.no.estex])
  dimnames(or.vcov) <- list(names(fixef(model)[original.no.estex]), 
                            names(fixef(model)[original.no.estex]))
  or.test <- coef(summary(model))[original.no.estex, 3]
  if (!obs) {
    if (is.null(select)) {
      alt.fixed <- matrix(ncol = n.pred, nrow = n.groups, 
                          data = NA)
      dimnames(alt.fixed) <- list(grouping.names, names(fixef(model))[original.no.estex])
      alt.se <- matrix(ncol = n.pred, nrow = n.groups, 
                       data = NA)
      dimnames(alt.se) <- list(grouping.names, names(fixef(model))[original.no.estex])
      alt.vcov <- list()
      alt.test <- matrix(ncol = n.pred, nrow = n.groups, 
                         data = NA)
      dimnames(alt.test) <- list(grouping.names, names(fixef(model))[original.no.estex])
      for (i in 1:n.groups) {
        if (count == TRUE) {
          print(n.groups + 1 - i)
        }
        model.updated <- exclude.influence(model = model, 
                                           grouping = group, level = grouping.names[i], 
                                           gf = gf, delete = delete)
        altered.no.estex <- which(substr(names(fixef(model.updated)), 
                                         1, 6) != "estex.")
        alt.fixed[i, ] <- as.matrix(fixef(model.updated)[altered.no.estex])
        alt.se[i, ] <- as.matrix(se.fixef(model.updated)[altered.no.estex])
        alt.vcov[[i]] <- as.matrix(vcov(model.updated)[altered.no.estex, 
                                                       altered.no.estex])
        alt.test[i, ] <- as.matrix(coef(summary(model.updated))[, 
                                                                3][altered.no.estex])
      }
    }
    if (!is.null(select)) {
      model.updated <- exclude.influence(model, group, 
                                         select, gf = gf, delete = delete)
      altered.no.estex <- which(substr(names(fixef(model.updated)), 
                                       1, 6) != "estex.")
      alt.fixed <- matrix(ncol = n.pred, nrow = 1, data = fixef(model.updated)[altered.no.estex])
      dimnames(alt.fixed) <- list("Altered model", 
                                  names(fixef(model.updated))[altered.no.estex])
      alt.se <- matrix(ncol = n.pred, nrow = 1, data = se.fixef(model.updated)[altered.no.estex])
      dimnames(alt.se) <- list("Altered model", names(fixef(model.updated))[altered.no.estex])
      alt.vcov <- list()
      alt.vcov[[1]] <- as.matrix(vcov(model.updated)[altered.no.estex, 
                                                     altered.no.estex])
      dimnames(alt.vcov[[1]]) <- list(names(fixef(model.updated)[altered.no.estex]), 
                                      names(fixef(model.updated)[altered.no.estex]))
      alt.test <- matrix(ncol = n.pred, nrow = 1, data = coef(summary(model.updated))[, 
                                                                                      3][altered.no.estex])
      dimnames(alt.test) <- list("Altered model", 
                                 names(fixef(model.updated))[altered.no.estex])
    }
  }
  if (obs) {
    if (is.null(select)) {
      alt.fixed <- matrix(ncol = n.pred, nrow = n.obs, 
                          data = NA)
      dimnames(alt.fixed) <- list(1:n.obs, names(fixef(model))[original.no.estex])
      alt.se <- matrix(ncol = n.pred, nrow = n.obs, data = NA)
      dimnames(alt.se) <- list(1:n.obs, names(fixef(model))[original.no.estex])
      alt.vcov <- list()
      alt.test <- matrix(ncol = n.pred, nrow = n.obs, data = NA)
      dimnames(alt.test) <- list(1:n.obs, names(fixef(model))[original.no.estex])
      for (i in 1:n.obs) {
        if (count == TRUE) {
          print(n.obs + 1 - i)
        }
        model.updated <- exclude.influence(model, obs = i)
        altered.no.estex <- which(substr(names(fixef(model.updated)), 
                                         1, 6) != "estex.")
        alt.fixed[i, ] <- as.matrix(fixef(model.updated)[altered.no.estex])
        alt.se[i, ] <- as.matrix(se.fixef(model.updated)[altered.no.estex])
        alt.vcov[[i]] <- as.matrix(vcov(model.updated)[altered.no.estex, 
                                                       altered.no.estex])
        alt.test[i, ] <- as.matrix(coef(summary(model.updated))[, 
                                                                3][altered.no.estex])
      }
    }
    if (!is.null(select)) {
      model.updated <- exclude.influence(model, obs = select)
      altered.no.estex <- which(substr(names(fixef(model.updated)), 
                                       1, 6) != "estex.")
      alt.fixed <- matrix(ncol = n.pred, nrow = 1, data = fixef(model.updated)[altered.no.estex])
      dimnames(alt.fixed) <- list("Altered model", 
                                  names(fixef(model.updated))[altered.no.estex])
      alt.se <- matrix(ncol = n.pred, nrow = 1, data = se.fixef(model.updated)[altered.no.estex])
      dimnames(alt.se) <- list("Altered model", names(fixef(model.updated))[altered.no.estex])
      alt.vcov <- list()
      alt.vcov[[1]] <- as.matrix(vcov(model.updated)[altered.no.estex, 
                                                     altered.no.estex])
      dimnames(alt.vcov[[1]]) <- list(names(fixef(model.updated)[altered.no.estex]), 
                                      names(fixef(model.updated)[altered.no.estex]))
      alt.test <- matrix(ncol = n.pred, nrow = 1, data = coef(summary(model.updated))[, 
                                                                                      3][altered.no.estex])
      dimnames(alt.test) <- list("Altered model", 
                                 names(fixef(model.updated))[altered.no.estex])
    }
  }
  estex <- list(or.fixed = or.fixed, or.se = or.se, or.vcov = or.vcov, 
                or.test = or.test, alt.fixed = alt.fixed, alt.se = alt.se, 
                alt.vcov = alt.vcov, alt.test = alt.test)
  class(estex) <- "estex"
  return(estex)
}
