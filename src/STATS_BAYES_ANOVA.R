#/***********************************************************************
# * Licensed Materials - Property of IBM 
# *
# * IBM SPSS Products: Statistics Common
# *
# * (C) Copyright IBM Corp. 2015
# *
# * US Government Users Restricted Rights - Use, duplication or disclosure
# * restricted by GSA ADP Schedule Contract with IBM Corp. 
# ************************************************************************/

# author__ = "IBM SPSS, JKP"
# version__ = "1.0.0"

# History
# 18-sep-2015 Original Version


gtxt <- function(...) {
    return(gettext(...,domain="STATS_BAYES_ANOVA"))
}

gtxtf <- function(...) {
    return(gettextf(...,domain="STATS_BAYES_ANOVA"))
}

kwdmap = list("allmodels"="all", "stepdown"="top", "stepup"="bottom", "allwithmain"="withmain")
### MAIN ROUTINE ###
doBayesanova = function(dep, indepfixed=NULL, indeprandom=NULL, models="allwithmain", 
    comparison=NULL, maxmodels=10000,
    plotbayesf=FALSE, index=NULL, rscalecontfixed="medium", rscalecontrandom="medium",
    omitposteriorrandom=TRUE, iterations=10000, bayesfactoriterations=10000,
    modelsource="none", modelfile=NULL, workspaceaction="clear", modelfileout=NULL) {
    # Estimate Bayes regression
    
    # The modelsource and modelfile
    # parameters are not implemented, awaiting requests for that functionality

    setuplocalization("STATS_BAYES_ANOVA")
    
    # A warnings proc name is associated with the regular output
    # (and the same omsid), because warnings/errors may appear in
    # a separate procedure block following the regular output
    procname=gtxt("Bayesian ANOVA")
    warningsprocname = gtxt("Bayesian ANOVA: Warnings")
    omsid="STATSBAYESANOVA"
    warns = Warn(procname=warningsprocname,omsid=omsid)

    tryCatch(library(BayesFactor), error=function(e){
        warns$warn(gtxtf("The R %s package is required but could not be loaded.", "cmprsk"),dostop=TRUE)
        }
    )
    if (!is.null(spssdictionary.GetWeightVariable())) {
        warns$warn(
            gtxt("The dataset is weighted, but case weights are not used in this procedure except for screening out cases with a nonpositive weight"),
            dostop=FALSE)
    }
    if (!is.null(spssdata.GetSplitVariableNames())) {
        warns$warn(
            gtxt("Split variables are not honored by this procedure"),
            dostop=FALSE)
    }
    if (is.null(c(indepfixed, indeprandom))) {
        warns$warn(gtxt("At least one independent variable must be specified"),
            dostop=TRUE)
    }
    if (length(intersect(indepfixed, indeprandom)) > 0) {
        warns$warn(gtxt("The same variable cannot be both a fixed and a random factor"),
            dostop=TRUE)
    }
    if (!is.null(comparison) && comparison == 0) {
        comparison = NULL
    }
    # Allow for estimating a single equation
    if (models == "single") {
        comparison = NULL
        index = 1
        plotbayesf = FALSE
    }
    alldata = c(dep, indepfixed, indeprandom)
    frml = paste(dep, paste(c(indepfixed, indeprandom), collapse="+", sep="+"), sep="~")
    rscalecontfixed = scales(rscalecontfixed, warns)
    rscalecontrandom = scales(rscalecontrandom, warns)
    allargs = as.list(environment())
    dta = spssdata.GetDataFromSPSS(alldata, missingValueToNA=TRUE,
        factorMode="levels")
    if (!all(as.logical(lapply(dta[-1],is.factor)))) {
        warns$warn(gtxt("Continuous variables cannot be used in this procedure"),
            dostop=TRUE)
    }
    if (is.factor(dta[[1]])) {
        warns$warn(gtxt("The dependent variable must have a continuous (scale) measurement level"),
            dostop=TRUE)
    }
    # The procedure does not allow missing values
    allargs$ncases = nrow(dta)
    dta = dta[complete.cases(dta),]
    allargs$nvalid = nrow(dta)

    if (models != "single") {
        res = tryCatch(anovaBF(as.formula(frml), whichRandom=indeprandom,
                data=dta, whichModels=kwdmap[models], iterations=bayesfactoriterations,
                progress=FALSE, rscaleFixed=rscalecontfixed, rscaleRandom=rscalecontrandom),
            error = function(e) {
                warns$warn(e$message, dostop=TRUE)
            }
        )
    } else {
        res = tryCatch(lmBF(as.formula(frml), whichRandom=indeprandom,
            data=dta, progress=FALSE, iterations=bayesfactoriterations,
            rscaleFixed=rscalecontfixed, rscaleRandom=rscalecontrandom),
            error = function(e) {
                warns$warn(e$message, dostop=TRUE)
            }
        )
    }
    if (!is.null(allargs$comparison)) {
        allargs$comparison = checkcomparison(allargs$comparison, res, warns)

        res = tryCatch(res/res[allargs$comparison],
            error = function(e) {warns$warn(e, dostop=TRUE)}
        )
    }

    post = doposterior(allargs, res, 
            omitposteriorrandom=omitposteriorrandom, indeprandom, warns)
    displayresults(allargs, res, post, warns)
    
    if (!is.null(modelfile)) {
        save(allargs, res, post, file=modelfile)
    }
    if (workspaceaction == "retain" && is.null(modelfile)) {
        assign("allargs", allargs, envir=.GlobalEnv)
        assign("res", res, envir=.GlobalEnv)
        assign("post", post, envir=.GlobalEnv)
    }
    warns$display()
}

checkcomparison = function(comparison, res, warns) {
    # check comparison spec and return if okay
    if (is.null(comparison)) {
        return(NULL)
    }
    if (comparison > length(res)) {
        warns$warn(gtxtf(
            "The comparison or index model number is greater than the number of models, which is %s. Substituting last model", 
            length(res)), dostop=FALSE)
        return(length(res))
    } else {
        return(comparison)
    }
}

slist = list("medium" = .7071, 'wide'=1.0, "ultrawide"=1.4142)
scales = function(scale, warns) {
    # return a numeric scale item, resolving certain strings

    if (is.null(scale)) {
        return(.7071)
    }
    scale = tolower(scale)
    if (!is.numeric(scale)) {
        scale = slist[[tolower(scale)]]
    }

    if (is.null(scale) || scale <= 0 ) {
        warns$warn(gtxt("An invalid value was given for the prior scale"),
                   dostop=TRUE)
    }
    return(scale)
}
doposterior = function(allargs, res, omitposteriorrandom, indeprandom, warns) {
    # calculate posterior distribution if model index specified
    # if omitposteriorrandom is TRUE, randompredictors are excluded

    if (is.null(allargs$index)) {
        return(NULL)
    }
    allargs$index = checkcomparison(allargs$index, res, warns)
    arglist = list(model=res, index=allargs$index, iterations=allargs$iterations,
        progress=FALSE)

    if (is.null(indeprandom) || !omitposteriorrandom) {
        filter = NULL
    } else {
        # build an re that will excluded all random factor values
        # (yes, this will rewrite the re multiple times, but that's not an issue here)
        # The filter has to match the variable name exactly without the category value.
        filter = ""
        for (v in indeprandom) {
            filter = paste(filter, sprintf("^%s$", c(v)), sep="|",collapse="")
        }
        arglist["columnFilter"]=substr(filter, 2, nchar(filter))
    }

    post = tryCatch(do.call(posterior, arglist),
                error=function(e) {warns$warn(e, dostop=TRUE)},
                warning = function(w) {warns$warn(e, dostop=FALSE)
                   return(NULL)}
    )
#     post = tryCatch(posterior(model=res, index=allargs$index, iterations=allargs$iterations,
#         , progress=FALSE),
#         error=function(e) {warns$warn(e, dostop=TRUE)},
#         warning = function(w) {warns$warn(e, dostop=FALSE)
#             return(NULL)}
#     )
    return(post)
}
    
scaletrans=list("medium"=gtxt("medium"), "wide"=gtxt("wide"), "ultrawide"=gtxt("ultrawide"))
waction=list("clear"="clear", "retain"="retain")

displayresults = function(allargs, res, post, warns) {
    # display results
    # allargs is the parameter set
    
    ressum = extractBF(res)

    StartProcedure(allargs[["procname"]], allargs[["omsid"]])
    
    # summary results
    # input specifications
    # although groups can be specified (cengroup), separate results are not
    # produced.
    lbls = c(gtxt("Dependent Variable"),
             gtxt("Fixed Factors"),
             gtxt("Random Factors"),
             gtxt("Comparison Model"),
             gtxt("Number of Cases"),
             gtxt("Number of Valid Cases"),
             gtxt("Prior Scale (Fixed"),
             gtxt("Prior Scale (Random)"),
             gtxt("Posterior Model Index"),
             gtxt("Posterior Iterations"),
             gtxt("Bayes Factor Integration Iterations"),
             gtxt("Workspace Action"),
             gtxt("Output Model File")
    )

    vals = c(
            allargs$dep,
            ifelse(is.null(allargs$indepfixed), gtxt("--NA--"), paste(allargs$indepfixed, collapse=" ")),
            ifelse(is.null(allargs$indeprandom), gtxt("--NA--"), paste(allargs$indeprandom, collapse=" ")),
            ifelse(is.null(allargs$comparison), 
                paste(gtxt("Intercept"), paste(allargs$indeprandom, collapse="+", sep="+"), collapse="+"), row.names(ressum)[allargs$comparison]),
            allargs$ncases,
            allargs$nvalid,
            allargs$rscalecontfixed,
            allargs$rscalecontrandom,
            ifelse(is.null(allargs$index), gtxt("--NA--"), allargs$index),
            ifelse(is.null(allargs$index), gtxt("--NA--"), allargs$iterations),
            allargs$bayesfactoriterations,
            waction[allargs$workspaceaction],
            ifelse(is.null(allargs$modelfile), gtxt("--NA--"), allargs$modelfile)
    )

    spsspivottable.Display(data.frame(cbind(vals), row.names=lbls), title = gtxt("Summary"),
        collabels=c(gtxt("Summary")), templateName="BAYESANOVASUMMARY", outline=gtxt("Bayes ANOVA Summary"),
        caption = gtxtf("Computations done by R package BayesFactor, version: %s", packageVersion("BayesFactor"))
    )

    bf = data.frame(seq(1: length(res)),ressum[1:2])
    bf[3] = bf[3] * 100.
    bf = data.frame(bf, length(res) - rank(bf[2]) + 1)
    # add in posterior probabilities excluding Intercept only model
    
    # construct posterior probabilities and merge with Bayes factors
    # The order for probabilities may not be the same as for the Bayes factors
    # which requires a few extra steps to get things merged
    # the BF data frame may not have the intercept row, so that row may be discarded
    postprob = data.frame(as.BFprobability(newPriorOdds(res) * res))[1]

    bf = merge(bf, postprob, by="row.names")
    bf = bf[order(bf[[2]]),]
    row.names(bf) = bf[["Row.names"]]
    bf = bf[-1]

    names(bf) = c(gtxt("Model Number"),
        gtxt("Bayes Factor"), gtxt("Error (+-%)"), gtxt("Rank"),
        gtxt("Posterior Probabilities (Equal Prior)"))

    if (allargs$models == "allmodels") {
        caption = gtxt("All models")
    } else if (allargs$models == "allwithmain") {
        caption = gtxt("All models that include corresponding main effects")
    } else if (allargs$models == "stepdown") {
        caption = gtxt("One variable at a time removed from full model")
    } else if (allargs$models == "stepup") {
        caption = gtxt("All single variable models")
    } else {
        caption = "Single model"
    }
    spsspivottable.Display(bf,
        title=gtxt("Bayes Factors"),
        rowdim=gtxt("Equation"), 
        hiderowdimtitle=FALSE,
        templateName="BAYESREGRFACTORS",
        outline=gtxt("Bayes Factors"),
        caption=caption
    )
    

    if (!is.null(allargs$index)) {
        postsum = summary(post)
        postsumstats = postsum$statistics[,-4]  # omit time series SEs
        # extras in rows won't conflict with variable names, because a category
        # suffix is always added, except that g_ variables are a little risky

        names(postsumstats) = c(gtxt("Mean"), gtxt("Std. Deviation"), gtxt("SE Mean"))
        spsspivottable.Display(
            postsumstats, 
            title=gtxtf("Posterior Summary Statistics for Model %s", allargs$index),
            rowdim=gtxt("Variables"),
            hiderowdimtitle=FALSE,
            templateName="BAYESREGRPOSTSTATS",
            outline=gtxt("Posterior Summary Statistics")
        )
        
        postsumquant = postsum$quantiles
        spsspivottable.Display(
            postsumquant,
            title=gtxtf("Posterior Quantiles for Model %s", allargs$index),
            rowdim=gtxt("Variables"),
            hiderowdimtitle=FALSE,
            templateName="BAYESREGRPOSTQUANTILES",
            outline=gtxt("Posterior Quantiles")
        )
    }

    if (allargs$plotbayesf) {
        plot(res)
    }
    
    spsspkg.EndProcedure()
}

Warn = function(procname, omsid) {
    # constructor (sort of) for message management
    lcl = list(
        procname=procname,
        omsid=omsid,
        msglist = list(),  # accumulate messages
        msgnum = 0
    )
    # This line is the key to this approach
    lcl = mylist2env(lcl) # makes this list into an environment

    lcl$warn = function(msg=NULL, dostop=FALSE, inproc=FALSE) {
        # Accumulate messages and, if dostop or no message, display all
        # messages and end procedure state
        # If dostop, issue a stop.

        if (!is.null(msg)) { # accumulate message
            assign("msgnum", lcl$msgnum + 1, envir=lcl)
            # There seems to be no way to update an object, only replace it
            m = lcl$msglist
            m[[lcl$msgnum]] = msg
            assign("msglist", m, envir=lcl)
        } 

        if (is.null(msg) || dostop) {
            lcl$display(inproc)  # display messages and end procedure state
            if (dostop) {
                stop(gtxt("End of procedure"), call.=FALSE)  # may result in dangling error text
            }
        }
    }
    
    lcl$display = function(inproc=FALSE) {
        # display any accumulated messages as a warnings table or as prints
        # and end procedure state, if any

        if (lcl$msgnum == 0) {   # nothing to display
            if (inproc) {
                spsspkg.EndProcedure()
            }
        } else {
            if (!inproc) {
                procok =tryCatch({
                    StartProcedure(lcl$procname, lcl$omsid)
                    TRUE
                    },
                    error = function(e) {
                        FALSE
                    }
                )
            }
            if (procok) {  # build and display a Warnings table if we can
                table = spss.BasePivotTable("Warnings ","Warnings") # do not translate this
                rowdim = BasePivotTable.Append(table,Dimension.Place.row, 
                    gtxt("Message Number"), hideName = FALSE,hideLabels = FALSE)

                for (i in 1:lcl$msgnum) {
                    rowcategory = spss.CellText.String(as.character(i))
                    BasePivotTable.SetCategories(table,rowdim,rowcategory)
                    BasePivotTable.SetCellValue(table,rowcategory, 
                        spss.CellText.String(lcl$msglist[[i]]))
                }
                spsspkg.EndProcedure()   # implies display
            } else { # can't produce a table
                for (i in 1:lcl$msgnum) {
                    print(lcl$msglist[[i]])
                }
            }
        }
    }
    return(lcl)
}

mylist2env = function(alist) {
    env = new.env()
    lnames = names(alist)
    for (i in 1:length(alist)) {
        assign(lnames[[i]],value = alist[[i]], envir=env)
    }
    return(env)
}

# localization initialization
setuplocalization = function(domain) {
    # find and bind translation file names
    # domain is the root name of the extension command .R file, e.g., "SPSSINC_BREUSCH_PAGAN"
    # This would be bound to root location/SPSSINC_BREUSCH_PAGAN/lang

    fpath = Find(file.exists, file.path(.libPaths(), paste(domain, ".R", sep="")))
    bindtextdomain(domain, file.path(dirname(fpath), domain, "lang"))
} 
# override for api to account for extra parameter in V19 and beyond
StartProcedure <- function(procname, omsid) {
    if (substr(spsspkg.GetSPSSVersion(),1, 2) >= 19) {
        spsspkg.StartProcedure(procname, omsid)
    }
    else {
        spsspkg.StartProcedure(omsid)
    }
}



Run = function(args) {
    #Execute the STATS COMPRISK command

    cmdname = args[[1]]
    args = args[[2]]
    oobj = spsspkg.Syntax(list(
        spsspkg.Template("DEPENDENT", subc="", ktype="existingvarlist", var="dep"),
        spsspkg.Template("INDEPFIXED", subc="", ktype="existingvarlist", var="indepfixed", islist=TRUE),
        spsspkg.Template("INDEPRANDOM", subc="", ktype="existingvarlist", var="indeprandom", islist=TRUE),
        spsspkg.Template("MODELS", subc="", ktype="str", var="models",
            vallist=list("allmodels", "allwithmain", "stepdown", "stepup", "single")),
        
        spsspkg.Template("COMPARISON", subc="OPTIONS", ktype="int", var="comparison"),
        spsspkg.Template("MAXMODELS", subc="OPTIONS", ktype="str", var="maxmodels"),
        spsspkg.Template("PLOTMODELS", subc="OPTIONS", ktype="bool", var="plotbayesf"),
        spsspkg.Template("POSTERIORINDEX", subc="OPTIONS", ktype="int", var="index"),
        spsspkg.Template("OMITPOSTERIORRANDOM", subc="OPTIONS", ktype="bool", var="omitposteriorrandom"),
        spsspkg.Template("ITERATIONS", subc="OPTIONS", ktype="int", var="iterations",
            vallist=list(2)),
        spsspkg.Template('BAYESFACTORITERATIONS', subc="OPTIONS", ktype="int", var="bayesfactoriterations",
            vallist=list(1)),
        spsspkg.Template("PRIORSCALEFIXED", subc="OPTIONS", ktype="literal", var="rscalecontfixed"),
        spsspkg.Template("PRIORSCALERANDOM", subc="OPTIONS", ktype="literal", var="rscalecontrandom"),
        
        spsspkg.Template("WORKSPACE", subc="SAVE", ktype="str", var="workspaceaction",
            vallist=list("retain", "clear")),
        spsspkg.Template("MODELFILE", subc="SAVE", ktype="literal", var="modelfileout")
    ))

    # A HELP subcommand overrides all else
    if ("HELP" %in% attr(args,"names")) {
        helper(cmdname)
    }
    else {
        res <- spsspkg.processcmd(oobj, args, "doBayesanova")
    }
}

helper = function(cmdname) {
    # find the html help file and display in the default browser
    # cmdname may have blanks that need to be converted to _ to match the file
    
    fn = gsub(" ", "_", cmdname, fixed=TRUE)
    thefile = Find(file.exists, file.path(.libPaths(), fn, "markdown.html"))
    if (is.null(thefile)) {
        print("Help file not found")
    } else {
        browseURL(paste("file://", thefile, sep=""))
    }
}
if (exists("spsspkg.helper")) {
assign("helper", spsspkg.helper)
}
