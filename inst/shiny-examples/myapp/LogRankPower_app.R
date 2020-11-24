#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(waiter)
library(shiny)
library(tidyverse)
library(DT)
library(lubridate)
#library(rio)
library(janitor)
library(shinyjs)
library(ggplot2)
library(shinythemes)
library(survival)
#library(Biobase)
library(rmarkdown)

# Define UI for application that draws a histogram
ui <- fluidPage(
  #use_waiter(), # include dependencies
    # Application title
    titlePanel("Log Rank Two Sample Power"),

    # Inputs
    sidebarLayout(
        sidebarPanel(
            numericInput("medH0",
                         "Median Survival time (months) in Control Group:",
                         min = 1,
                         max = 120,
                         step = 1,
                         value = 5),

            numericInput("medH1",
                         "Median Survival time (months) in Treatment Group:",
                         min = 1,
                         max = 120,
                         step = 1,
                         value = 9),

            numericInput("n.sim",
                        "number of simulations:",
                        min = 1,
                        max = 10000,
                        step = 100,
                        value = 1000),

            numericInput("n0",
                         "number of patients in the control group:",
                         min = 1,
                         max = 1000,
                         step = 24,
                         value = 50),

            numericInput("n1",
                         "number of patients in the treatment group:",
                         min = 1,
                         max = 1000,
                         step = 25,
                         value = 50),


            numericInput("time.accrual",
                         "Accrual time (months):",
                         min = 1,
                         max = 120,
                         step = 1,
                         value = 24),

            numericInput("time.followup",
                         "Follow up time (months):",
                         min = 1,
                         max = 120,
                         step = 1,
                         value = 36),

            radioButtons("con.CI", "Confidence interval:",
                          choiceNames = list("95%" , "90%"  ),
                         choiceValues  = list( 0.95  , 0.90  ) ),
            # Show data table --------------------------------------------------------
            checkboxInput(inputId = "showdata",
                          label = "Show data table",
                          value = FALSE),

            actionButton("goButton", "Calculate!" ) ,

            radioButtons('format', 'Document format', c('Word', 'PDF'),

                         inline = TRUE),

            downloadButton('downloadReport')

        ), #end o'sidebarPanel

        # Show a plot of the generated distribution
        mainPanel(
          tags$style(type='text/css', '#overone {color: white;}'),

          h1("Summary"),
          textOutput("summary"),
          conditionalPanel( condition = 'output.overone == "true"',
                            textOutput("testingoveronesummary")  ),
          conditionalPanel( condition = 'output.overone != "true"',
                            textOutput("testingoveronesummaryfalse")  ),
          h2("Study design"),
          textOutput("studydesign"),

          h2("Statistical power"),
          h3("Power and Type I error"),
          textOutput("power"),

          h3("Distribution of HR"),
          textOutput("Relationship"),
          textOutput("HRdistribution"),
          plotOutput("boxplots"),

          h3("Relationship of HR and p value"),
          textOutput("relationship"),
          plotOutput("scatterplot"),

          plotOutput("scatterplot2"),

          h3("Relationship of median survival time between treatment and control stratified by HR"),

          textOutput("relationshipratio"),
          plotOutput("scatterplot3"),
          br(), br(),

          textOutput("overone" ),

          conditionalPanel( condition = 'output.overone == "true"',
                            textOutput("testingoverone"),
                            plotOutput("kmplot")
                            ),


          DT::dataTableOutput(outputId = "outputtable")


        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  #w <- Waiter$new()

      outcomesdata <- eventReactive(input$goButton, {
        input$goButton
        #w$show()
     out <- function(medH0 =  input$medH0,
                 medH1 =  input$medH1,
                 n0 =  input$n0,
                 n1 =  input$n1,
                 time.accrual =  input$time.accrual,
                 time.followup =  input$time.followup,
                 n.sim =  input$n.sim,
                 con.CI =  as.numeric(input$con.CI)
        ) {
            #  require(survminer)
            hazard.rate0 <- log(2) / medH0
            med0 <- median(rexp(n.sim, rate = hazard.rate0))
            hazard.rate1 <- log(2) / medH1
            med1 <- median(rexp(n.sim, rate = hazard.rate1))

            time.total <- time.accrual + time.followup
            event.rate.H0 <- 1 - 1 / ((time.followup) * hazard.rate0) * (exp(-time.accrual * hazard.rate0) - exp(-time.total * hazard.rate0))
            event.rate.H1 <- 1 - 1 / ((time.followup) * hazard.rate1) * (exp(-time.accrual * hazard.rate1) - exp(-time.total * hazard.rate1))

            HRatio <- hazard.rate1/hazard.rate0
            confirmHR<-medH0/medH1
              # print(c(HRatio = HRatio, confirmHR = confirmHR , HR0=hazard.rate0, med0=med0, er0=event.rate.H0,
              #            HR1=hazard.rate1, med1=med1, er1=event.rate.H1,
              #            totaltime = time.total, acctime = time.accrual, futime = time.followup))


            p.tmp <- matrix(numeric(n.sim * 20), n.sim, 20)
            data.list.tmp <- as.list(numeric(n.sim))
            for (i in 1:n.sim)
            { #print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ loop ")
                x0 <- rexp(n0, rate = hazard.rate0)
                c0 <- runif(n0, time.accrual, time.total)
                x1 <- rexp(n1, rate = hazard.rate1)
                c1 <- runif(n1, time.accrual, time.total)
                x10 <- c(x0, x1)
                c10 <- c(c0, c1)
                data1 <- data.frame(S_time = x10, C_time = c10, Obs_time = ifelse(x10 > c10, c10, x10),
                                    censor = ifelse(x10 > c10, 0, 1),
                                    group = rep(c("C", "T"), c(n0, n1)))
                # L1=(sum(data1$censor)-event.rate.H0*n)/sqrt(event.rate.H0*n)
                s1 <- Surv(data1$Obs_time, data1$censor)
                fit1 <- survfit(s1 ~ group, conf.int = con.CI, data = data1)
                cox1 <- coxph(s1 ~ group, data = data1)
                #print(summary(cox1))
                # surv.med<-summary(fit1)$table[c("median" ,    "0.9LCL" ,    "0.9UCL")]
                surv.med <- summary(fit1)$table[, -(1:6)]
                surv.med <- c(surv.med[1, ], surv.med[2, ]) # MED and CI from KM fit
                names(surv.med) <- paste(rep(c("C", "T"), each = 3), names(surv.med), sep = "_")
                coef1 <- coef(summary(cox1))[1, ]
                #print(names(summary(cox1)))
                #print(summary(cox1))
                test.stat <- unlist(summary(cox1)[c("sctest", "logtest", "waldtest")])
                #print(unlist(summary(cox1)))

                tmp1 <- c(surv.med, coef1, test.stat)
                p.tmp[i, ] <- tmp1
                #print(p.tmp)
                data.list.tmp[[i]] <- data1
            }
            colnames(p.tmp) <- names(tmp1)
            p.tmp <- data.frame(p.tmp) %>% mutate(ratiomed = C_median/T_median,
                                difference = C_median-T_median,
                                HRnum =  exp.coef.)
            list(data = data.list.tmp, ans = data.frame(p.tmp))
            #return(data.frame(p.tmp))
     }
     #Sys.sleep(3)
     #w$hide()
     out()


    })

      output$outcomesdata <- renderTable({
          outcomesdata()[[2]]

      })

      output$outcomesrawdata <- renderTable({
        outcomesdata()[[1]][[1]]

      })

      output$outputtable <- DT::renderDataTable({
        if(input$showdata){
          DT::datatable(data = outcomesdata()[[2]],
                        options = list(pageLength = 10),
                        rownames = FALSE)
        }
      })

      a2 <- reactive({outcomesdata()[[2]] %>%
            #filter((outcomesdata()[[2]]$sctest.pvalue < 0.05) & (outcomesdata()[[2]]$coef < 0)) %>%
            filter((outcomesdata()[[2]]$sctest.pvalue < 0.05) & (outcomesdata()[[2]]$coef < 0)) %>%
            mutate(HR = cut(exp.coef., breaks = c(floor(min(exp.coef.)*10):ceiling(max(exp.coef.)*10)) / 10))
        })

      output$a2 <- renderTable({
        a2()
       })

       categories <- reactive({ names(table(a2()$HR))    })
       counts <- reactive({ unname(round(table(a2()$HR)/sum(table(a2()$HR))*100, 2))       })
       percents <- reactive({ unname(round(table(a2()$HR)/sum(table(a2()$HR))*100, 2))        })
       statement <- reactive({  paste(categories()," ",percents(), "%",sep='', collapse = ', ')         })
       POW <- reactive({ sum(I(outcomesdata()[[2]]$sctest.pvalue<0.05))/input$n.sim       })
       numberoutliers <- reactive({   numberoutliers <- dim(outcomesdata()[[2]] %>% filter(ratiomed>1 & HRnum<1))[1]  })
       percentoutliers <- reactive({  dim(outcomesdata()[[2]] %>% filter(ratiomed>1 & HRnum<1))[1]/input$n.sim*100 })
       # output$overone <- reactive({  I(NROW(outcomesdata()[[2]] %>% filter(ratiomed>1)) >0)  })
       # overonetf <- reactive({  I(NROW(outcomesdata()[[2]] %>% filter(ratiomed>1)) >0)       })

         output$overone <- reactive({  I(NROW(outcomesdata()[[2]] %>% filter(sctest.pvalue<0.05 & ratiomed>1)) >0)  })
         overonetf <- reactive({  I(NROW(outcomesdata()[[2]] %>% filter(sctest.pvalue<0.05 & ratiomed>1)) >0)       })


         #observe({ print(names(outcomesdata()[[2]]))          })

       output$summary <- renderText({
         paste(paste("	Summary:"," ",sep = "\n"),
               paste(
                 "This is a two-arm randomized trial with survival as primary endpoint.
 The study includes an enrollment time of ", input$time.accrual,"
 months and a follow-up time of at least ", input$time.followup," months per patient.
 The primary hypothesis is that the treatment will yield an improved median
 survival time of ", input$medH1 - input$medH0, "months compared to the
 control group of a ", input$medH0,"-month median survival time.   The control group has ", input$n0 ,"
 patients and the treatment group has", input$n1 ," for a total of", input$n0 + input$n1, "patients.",
 "	The proposed design with n = ",input$n0," per group will achieve ",
 sum(I(outcomesdata()[[2]]$sctest.pvalue<0.05))/input$n.sim," power to detect a hazard ratio (HR)
 of ",   round(input$medH0/input$medH1,3) ," controlled at a one-sided ",
                 (1-as.numeric(input$con.CI))*100 ,"% type I error based on log-rank test.",
 "If the treatment is effective with an effect size of HR = ",round(input$medH0/input$medH1,3) ,"(i.e., alternative
 hypothesis is true),  the estimated median HR will be on the target,", round(input$medH0/input$medH1,3) ,",
 with one-sided 95% confidence interval (95% CI) of ",

   round(quantile(round(outcomesdata()[[2]]$HRnum,5), probs = c(0.025)),2),"-",
   round(quantile(round(outcomesdata()[[2]]$HRnum,5), probs = c(0.975)),2),

 "
 .The significant HR (log-rank p<0.05) ranges",
       round(min(outcomesdata()[[2]]$exp.coef.[outcomesdata()[[2]]$sctest.pvalue < 0.05]),2) ,
 "to", round(max(outcomesdata()[[2]]$exp.coef.[outcomesdata()[[2]]$sctest.pvalue < 0.05]),2) ,".
 HR and p value are highly correlated with low HR in small p value spectrum (r = ",
 round(cor((outcomesdata()[[2]]$coef), log(outcomesdata()[[2]]$sctest.pvalue)), 2) ," in logarithm scale). ",
 "Distribution of median survival time (95% CI) is ",
 round(quantile(round(outcomesdata()[[2]]$C_median,5), probs = c(0.025)),2),"-",
 round(quantile(round(outcomesdata()[[2]]$C_median,5), probs = c(0.975)),2)," months with median value of ", input$medH0,"
 months in the control group and ",
 round(quantile(round(outcomesdata()[[2]]$T_median,5), probs = c(0.025)),2),"-",
 round(quantile(round(outcomesdata()[[2]]$T_median,5), probs = c(0.975)),2)," months with median value of ", input$medH1," months in the treatment group.
 The 95% CI of median survival ratio is
 ",
 round(quantile(round(outcomesdata()[[2]]$C_median/outcomesdata()[[2]]$T_median,5), probs = c(0.025)),2),"-",
 round(quantile(round(outcomesdata()[[2]]$C_median/outcomesdata()[[2]]$T_median,5), probs = c(0.975)),2)," with median value
 of
 ",round(summary(outcomesdata()[[2]]$C_median/outcomesdata()[[2]]$T_median)[3],2),". For significant cases
 with log-rank p < 0.05, the range of median survival time of control and treatment groups and median survival ratio is",
 round(min(outcomesdata()[[2]]$C_median),2) ,
 "to", round(max(outcomesdata()[[2]]$C_median),2) ," months,",
 round(min(outcomesdata()[[2]]$T_median),2) ,
 "to", round(max(outcomesdata()[[2]]$T_median),2) ," months,",
 round(min(outcomesdata()[[2]]$C_median/outcomesdata()[[2]]$T_median),2) ,
 "to", round(max(outcomesdata()[[2]]$C_median/outcomesdata()[[2]]$T_median),2) ,", respectively.
 Median survival ratio is also positively correlated with
 HR, with a lower median survival ratio corresponds to a lower HR ( r = ",
 round(cor( outcomesdata()[[2]]$C_median/outcomesdata()[[2]]$T_median, outcomesdata()[[2]]$HRnum), 2),")."


               )
         )

       })

       output$testingoveronesummary <- renderText({

         paste("However, there are ",numberoutliers()," outliers (",percentoutliers(),
               "%) with opposite trend with a HR < 1 but \n with a median survival ratio > 1.
                          Close examination indicates most cases show \n treatment delay effect", sep = "")

       })

       output$testingoveronesummaryfalse <- renderText({
         if(numberoutliers()>0){
         paste("However, there are ",numberoutliers()," outliers (",percentoutliers(),
               "%) with opposite trend with a HR < 1 but \n with a median survival ratio > 1.
                          Yet none are significant at the alpha = ",1-as.numeric(input$con.CI)," level", sep = "")
         } else {paste(" ", sep = "")  }
       })


       output$studydesign <- renderText({
        paste(paste("	Study design:"," ",sep = "\n"),
        paste(
"This is a two-arm randomized trial with survival as primary endpoint.
The study includes an enrollment time of ", input$time.accrual,"
 months and a follow-up time of at least ", input$time.followup," months per patient.
 The primary hypothesis is that the treatment will yield an improved median
 survival time of ", input$medH1 - input$medH0, "months compared to the
 control group of a ", input$medH0,"-month median survival time.   The control group has ", input$n0 ,"
patients and the treatment group has", input$n1 ," for a total of", input$n0 + input$n1, "patients."
 )
)
      })

       output$power <- renderText({
        paste("	The proposed design with n = ",input$n0," per group will achieve ",
              sum(I(outcomesdata()[[2]]$sctest.pvalue<0.05))/input$n.sim," power to detect a hazard ratio (HR)
 of ",   round(input$medH0/input$medH1,3) ," controlled at a one-sided ",
              (1-as.numeric(input$con.CI))*100 ,"% type I error based on log-rank test.", sep = "")
      })

       output$HRdistribution <- renderText({
        paste("Simulation shows a median HR of ",   round(input$medH0/input$medH1,2) ," that ranges
from ", round(min(outcomesdata()[[2]]$exp.coef.), 2)," to ",round(max(outcomesdata()[[2]]$exp.coef.), 2),"
 (Figure 1a).  Distribution of significant HR (log-rank p<0.05)
 ranges ",  round(min(outcomesdata()[[2]]$exp.coef.[outcomesdata()[[2]]$sctest.pvalue < 0.05]),2)," to ",
               round(max(outcomesdata()[[2]]$exp.coef.[outcomesdata()[[2]]$sctest.pvalue < 0.05]),2) ," with ",
              round(summary(outcomesdata()[[2]]$exp.coef.[outcomesdata()[[2]]$sctest.pvalue < 0.05])[2],2),", ",
 round(summary(outcomesdata()[[2]]$exp.coef.[outcomesdata()[[2]]$sctest.pvalue < 0.05])[3],2), " and",
 round(summary(outcomesdata()[[2]]$exp.coef.[outcomesdata()[[2]]$sctest.pvalue < 0.05])[5],2) ,"
              in the 1st, 2nd, and 3rd quartiles, respectively (Figure 1b). ")
  })

       output$relationship <- renderText({
        paste("HR and p value show a high correlation (r = ",round(cor((outcomesdata()[[2]]$coef), log(outcomesdata()[[2]]$sctest.pvalue)), 2),";
 Figure 2) in logarithm scale with low HR in small p value spectrum")
      })

       output$relationshipratio <- renderText({
        paste("The statisically significant HR is grouped in the following increasing ordered categories: ",
                statement(),
              "; For HR < 0.5, the median survival time in the treatment group is always greater than
              the control group as shown in Figure 4 with all data points
              below the 45 degree red-color line (i.e., median survival ratio < 1). As HR exceeds 0.5,
              the number of intances in which the median survival ratio > 1 starts to increase."
              ,sep="" )
      })

       output$testing <- renderText({
        paste( POW(),numberoutliers(), percentoutliers(), overonetf() ,sep="  " )
       })

       output$testingoverone <- renderText({

           paste("However, there are ",numberoutliers()," outliers (",percentoutliers(),
                          "%) with opposite trend with a HR < 1 but \n with a median survival ratio > 1.
                          Close examination indicates most cases show \n treatment delay effect (Figure",
                          " 5",").", sep = "")

        })



      #---output report----

       output$downloadReport <- downloadHandler(

         filename = function() {

           paste('LogRankPowerReport', sep = '.',

                 switch(input$format, PDF = 'pdf', HTML = 'html', Word = 'docx'))

         },

         content = function(file) {

           out <- rmarkdown::render(input = 'LogRankPowerReport.Rmd',

                                    output_format =

                                      switch(input$format,

                                             PDF = rmarkdown::pdf_document(),

                                             HTML = rmarkdown::html_document(),

                                             Word = rmarkdown::word_document()

                                      ),

                                    params = list(set_title = input$project_title, set_author = input$author_input)

           )

           file.rename(out, file)

         }

       )


      # PLOTS ####
      output$scatterplot <- renderPlot({
        # w$show()
        plot((outcomesdata()[[2]]$coef), log(outcomesdata()[[2]]$sctest.pvalue), xlab = "log(HR)", ylab = "log(p)",
             main = paste("Fig 2. Correlation Coefficient Estimate =", round(cor((outcomesdata()[[2]]$coef), log(outcomesdata()[[2]]$sctest.pvalue)), 2)),
             cex.main = 1)
        abline(h = log(0.05), col = 2, lty = 2)
        axis(4, log(0.05), "p=0.05")
        # Sys.sleep(3)
        # w$hide()
    })

      output$boxplots <- renderPlot({


        if(FALSE){

         par(mfrow=c(1,2))
          boxplot(outcomesdata()[[2]]$exp.coef., ylab = "HR", main = "Fig 1a. Distribution of HR under H1", cex.main = .65,
                  ylim = c(min(outcomesdata()[[2]]$exp.coef.),max(outcomesdata()[[2]]$exp.coef.)))
          boxplot(outcomesdata()[[2]]$exp.coef.[outcomesdata()[[2]]$sctest.pvalue < 0.05], ylab = "HR (p < 0.05)",
                  main = "Fig 1b. Distribution of Significant HR under H1", cex.main = .65,
                  ylim = c(min(outcomesdata()[[2]]$exp.coef.),max(outcomesdata()[[2]]$exp.coef.)))
        } else {

          par(mfrow=c(1,2))
          boxplot(outcomesdata()[[2]]$exp.coef., ylab = "HR", main = "Fig 1a. Distribution of HR under H1", cex.main = .65,
                  ylim = c(min(outcomesdata()[[2]]$exp.coef.),max(outcomesdata()[[2]]$exp.coef.)))
          boxplot(outcomesdata()[[2]]$exp.coef.[outcomesdata()[[2]]$sctest.pvalue < 0.05], ylab = "HR (p < 0.05)",
                  main = "Fig 1b. Distribution of Significant HR under H1", cex.main = .65 )

        }



      })

      output$scatterplot2 <- renderPlot({

        plot(a2()$C_median / a2()$T_median, a2()$exp.coef., ylab = "Significant HR (p<0.05)",
             xlab = "Ratio of median survival time (control/treatment)",
             main = paste("Fig 3. Correlation Coefficient Estimate =", round(cor(a2()$C_median / a2()$T_median,
                                                                          log(a2()$sctest.pvalue)), 2)), cex.main = 1)
        abline(0, 1, col = 2)
        abline(v = 1, col = 3)
      })

      output$scatterplot3 <- renderPlot({
      a2() %>%
        ggplot(aes(x = T_median, y = C_median, colour = HR, shape = HR)   ) +
            labs(title = "Fig 4. Median survival time between treatment and control stratified by (significant) HR",
                 x = "Median of treatment group", y = "Median of control group") +
        geom_point() +
        theme(legend.position = "top") +
        #  geom_smooth(method="lm",colour='black')+
        geom_abline(slope = 1, intercept = 0, colour = "red") +
        #  geom_text(aes(min(T_median)*1.4, max(C_median), label = paste('Int=',intercept,"Slope =", slope, "\n")))+
        facet_wrap(~HR)})

      output$kmplot <- renderPlot({

        index1 <- (1:dim(outcomesdata()[[2]])[1])[(outcomesdata()[[2]]$sctest.pvalue < 0.05) & (outcomesdata()[[2]]$coef < 0)]
        index2 <- index1[(a2()$sctest.pvalue < 0.05) & (a2()$ratiomed > 1)  & (a2()$ratiomed == max(a2()$ratiomed))]

        data2 <- outcomesdata()[[1]][[index2[1]]]

        f1 <- survfit(Surv(data2$Obs_time, data2$censor) ~ data2$group)
        f2 <- coxph(Surv(data2$Obs_time, data2$censor) ~ data2$group)

        # DT::datatable(data = data2,
        #               options = list(pageLength = 10),
        #               rownames = FALSE)

        plot(f1, col = 1:2, lty = 1:2, xlab = "time", ylab = "survival probability")
        legend(10, 1, c("C", "T"), col = 1:2, lty = 1:2)
        med.tmp <- round(summary(f1)$table[, "median"], 2)
        title(paste("Fig 5. Estimated HR=", round(coef(summary(f2))[, "exp(coef)"], 2), "\n  Estimated median survival ratio=", round(med.tmp[1] / med.tmp[2], 2), "\n (", med.tmp[1], " in control versus ", med.tmp[2], " in treatment)", sep = ""),
              cex.main = 1)
        abline(h = 0.5, col = 3, lty = 2)
        axis(4, 0.5, 0.5)



      })

}

# Run the application
shinyApp(ui = ui, server = server)
