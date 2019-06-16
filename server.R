##---- SERVER ----
server <- function(input, output, session) {
    
    ## load data file  ----
    data_file <- reactive({
        req(input$file)
        data <- fread(input$file$datapath)
        return(data)
    })
    
    ## make data table and set 'head' function ---- 
    output$table <- renderTable({
        req(data_file)
        data <- data_file()
        if(input$dispout == "head") {
            return(head(data, 5))
        }
        else {
            return(data)
        }
    },caption = 'Raw Data Output Table', caption.placement = getOption("xtable.caption.placement", "top"))
    
    ## Head of raw data to visualise headers ----
    output$tablehead <- renderTable({
        req(data_file)
        data <- data_file()
        return(head(data, 3))
    },caption = 'Raw Data Output Table', caption.placement = getOption("xtable.caption.placement", "top"))
    
    ## Delta Ct Table ---- 
    output$dct_table <- renderTable({
        req(deltasummary)
        data <- deltasummary()
        colnames(data) <- c('Treatment','Gene','Delta Ct','Sample Size','Std Err','Transcript Abundance','Percent TrA','Lower TrA Error Pos','Upper TrA Error Pos')
        if(input$disp == "head") {
            return(head(data, 5))
        }
        else {
            return(data)
        }
    },caption = 'Delta Ct Summary', caption.placement = getOption("xtable.caption.placement", "top"))
    ## DeltaDelta Ct Table ----
    output$ddct_table <- renderTable({
        req(deltasummary)
        data <- ddCT_file()
        colnames(data) <- c('Treatment','Gene','Delta Delta Ct','Sample Size','Std Err','Fold Change','Fold Change Std Err','Lower FC Error Pos','Upper FC Error Pos')
        if(input$disp == "head") {
            return(head(data, 5))
        }
        else {
            return(data)
        }
    },caption = 'Delta Delta Ct Summary', caption.placement = getOption("xtable.caption.placement", "top"))
    
    
    # Fill the 'column' dropdown ----
    observeEvent( 
        input$file,{
            updateSelectInput(session,
                              'genecol',
                              choices=colnames(data_file()))
            updateSelectInput(session,
                              'treatcol',
                              choices=colnames(data_file()))
            updateSelectInput(session,
                              'ctcol',
                              choices=c('',colnames(data_file())))
            updateSelectInput(session,
                              'repcol',
                              choices=colnames(data_file()))
            output$selectedcol <- renderText({ 
                paste("Gene Column:", input$genecol)})
            output$selectedgene <- renderText({
                paste('Gene:', input$refgene)})
        })
    observeEvent(
        input$genecol,{
            updateSelectInput(session,
                              'refgene',
                              choices=c(unique(as.character(data_file()[[input$genecol]]))))
        })
    observeEvent(
        input$treatcol,{
            updateSelectInput(session,
                              'controlsample',
                              choices=c(unique(as.character(data_file()[[input$treatcol]]))))
        })
    
    ## delta CT data ----  
    deltasummary <- reactive({
        req(input$ctcol)
        df <- data_file()
        df <- df %>% rename('CT'=input$ctcol,
                            'Treatment'=input$treatcol,
                            'Gene'=input$genecol,
                            'Replicate'=input$repcol)
        df$CT<- as.numeric(as.character(df$CT))
        refname <- input$refgene
        df <- df %>% filter(!Task == 'NTC')
        df <- df %>% filter(!CT == 'Undetermined')
        
        df1 <- df %>%
            group_by(Replicate, Treatment, Gene) %>% 
            summarise(delta.CT = mean(CT)) %>% 
            mutate_at(vars(matches("delta")), funs(.- .[Gene == refname]))%>%
            filter(Gene != refname) %>% 
            group_by(Treatment, Gene) %>% 
            summarise_each(funs(mean(., na.rm = T), 
                                n = sum(!is.na(.)), 
                                se = sd(., na.rm = T)/sqrt(sum(!is.na(.)))), 
                           delta.CT) %>% 
            mutate(TrA = 2^-mean, # transcript abundance (TrA)
                   perc.TrA = TrA * 100, # percentage TrA
                   lw.TrA = 2^ - (mean + se),  # lower error position
                   up.TrA = 2^ - (mean - se)) %>%  # upper error position
            rename(dCT = mean) %>% 
            mutate_if(is.numeric, round, 5)
        df1
    })
    
    ##ddct data ----  
    ddCT_file <- reactive({
        req(input$ctcol)
        df <- deltasummary()
        refname <- input$refgene
        cont_samp <- input$controlsample
        
        df2 <- df %>% #####
        group_by(Gene) %>%
            mutate_at(vars(matches("dCT")), funs(.- .[Treatment == cont_samp])) %>% #specify control treatment
            select(1:5) %>% 
            rename(ddCT = dCT) %>% 
            mutate(fc = 2^-ddCT, # fold change ddCT
                   fc.se = se^2) %>% # SE fold change
            mutate_at(vars(matches("fc.se")), funs(sqrt(.+. [Treatment == cont_samp]))) %>% 
            mutate(lw.fc = 2^ - (ddCT + fc.se),  # lower error position
                   up.fc = 2^ - (ddCT - fc.se)) %>% # upper error position
            mutate_if(is.numeric, round, 5)
        df2
    })
    
    ##Create dCt Plot ----
    dctplot <- reactive({
        req(input$ctcol)
        data <- deltasummary()
        p <- ggplot(data = data, 
                    aes(x = Gene, 
                        y = TrA, 
                        fill = Treatment)) +
            geom_bar(stat = "identity", position = position_dodge(width = 0.7), 
                     colour = "black", width=0.6) +
            geom_errorbar(aes(ymin = lw.TrA, ymax = up.TrA, width = 0.1), # add standard error bars
                          position = position_dodge(width = 0.7)) +
            ylab(paste('Transcript Abundance relative to',input$refgene)) +
            ggtitle('Transcript Abundance plot') +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  legend.position="right",
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "black")) +
            scale_y_continuous(expand = c(0, 0))
    })
    #Output dCt Plot ----
    output$dctplot <- renderPlot({
        req(dctplot)
        dctplotPrint <- dctplot()
        dctplotPrint
    })
    ##Create ddCt Plot ---- 
    ddctplot <- reactive({
        req(input$ctcol)
        data <- ddCT_file()
        p <- ggplot(data = data, 
                    aes(x = Gene, 
                        y = fc, 
                        fill = Treatment)) +
            
            geom_bar(stat = "identity", position = position_dodge(width = 0.7), 
                     colour = "black", width=0.6) +
            
            geom_errorbar(aes(ymin = lw.fc, ymax = up.fc, width = 0.1), # add standard error bars
                          position = position_dodge(width = 0.7)) +
            ylab(paste('Fold Change relative to treatment',input$controlsample)) +
            ggtitle('Fold Change Plot') +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  legend.position="right",
                  panel.background = element_blank(),
                  axis.line = element_line(colour = "black")) +
            scale_y_continuous(expand = c(0, 0))
    })
    ##Output ddCt Plot ----
    output$ddctplot <- renderPlot({
        req(ddctplot)
        ddctplotPrint <- ddctplot()
        ddctplotPrint
    })
    ##Download DCt Plot ----  
    output$downloadPlot1 <- downloadHandler(
        filename = function(){paste('DeltaCtPlot.png')},
        content = function(file) {
            ggsave(file, plot = dctplot(), device = "png", width=8, height=6, units='in')
        }
    )
    ##Download DDCt Plot ----  
    output$downloadPlot2 <- downloadHandler(
        filename = function(){paste('DeltaDeltaCtPlot.png')},
        content = function(file) {
            ggsave(file, plot = ddctplot(), device = "png", width=8, height=6, units='in')
        }
    )
    ##Download DCt Table ----
    output$dct_table_DL <- downloadHandler(
        filename = function(){paste('DeltaCtSummary.csv')},
        content = function(file) {
            write.csv(deltasummary(), file, row.names = FALSE)
        }
    )
    ##Download DDCt Table ----
    output$ddct_table_DL <- downloadHandler(
        filename = function(){paste('DeltaDeltaCtSummary.csv')},
        content = function(file) {
            write.csv(ddCT_file(), file, row.names = FALSE)
        }
    )
    
} # end server
