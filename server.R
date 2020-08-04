shinyServer(function(input, output, session) {
    
    output$simulationUI<-renderUI({conditionalPanel(
        condition="input.simulation=='Yes'",
        fileInput('detectionThresholdsFile',label = "Upload detection thresholds file")
    )})
    
    output$simulationUI1<-renderUI({conditionalPanel(
        condition="input.simulation=='No'",
        radioButtons('exactSwitch',"How shall Quality Analysis be performed?",choices = c("Exact","Simulation"),selected = "Simulation",inline = T)
    )})
    
    output$simulationUI1.2<-renderUI({conditionalPanel(
        condition="input.simulation=='No'&&input.exactSwitch=='Exact'",
        h5("Note: selection of 'Exact' may lead to prolonged runtime."),
        hr(),
        h4("A: Evaluate Cell Fractions:"),
        radioButtons('simulate_germline_vaf',label="VAF for polymorphisms in the control sample",
                     choices = c("Exact (variant calling results)","Simulated (expected value 0.5)"),
                     selected = "Exact (variant calling results)",inline = T),
        sliderInput('simulate_rates',"Fractions to consider [%]",min=1,max=100,value=c(5,100)),
        numericInput('simulate_steps',"In steps of...",min=1,max=100,value=5),
        hr(),
        h4("B: Evaluate Window Sizes:"),
        sliderInput('window_percentile',"Percentile to consider for window selection [%]",min=0,max=100,value=95),
        hr(),
        radioButtons('simulation_optimization1',label="Select strategy for detection threshold optimization",
                     choices=c("Compromize (cell fraction and window size)","Force (cell fraction)","Force (window size)"),
                     selected= "Compromize (cell fraction and window size)",inline = F)
    )})
    
    output$simulationUI2<-renderUI({conditionalPanel(
        condition="input.simulation=='No'&&input.exactSwitch=='Simulation'",
        numericInput('simulations_nr',"Number of simulations to perform [10, 99999]",min=10,max=99999,value=500),
        hr(),
        h4("A: Evaluate Cell Fractions:"),
        sliderInput('simulate_rates',"Fractions to consider [%]",min=1,max=100,value=c(5,100)),
        numericInput('simulate_steps',"In steps of...",min=1,max=100,value=5),
        hr(),
        h4("B: Evaluate Window Sizes:"),
        sliderInput('window_percentile',"Percentile to consider for window selection [%]",min=0,max=100,value=95),
        hr(),
        radioButtons('simulation_optimization',label="Select strategy for detection threshold optimization",
                     choices=c("Compromize (cell fraction and window size)","Force (cell fraction)","Force (window size)"),
                     selected = "Compromize (cell fraction and window size)",inline = F)
    )})
    output$simulationUI3<-renderUI({conditionalPanel(
        condition="input.simulation_optimization=='Force (cell fraction)'||input.simulation_optimization1=='Force (cell fraction)'",
        numericInput('simulations_force_rate',"Minimum (possible) cell fraction to consider [%]",min=1,max=100,value=10)
    )})
    output$simulationUI4<-renderUI({conditionalPanel(
        condition="input.simulation_optimization=='Force (window size)'||input.simulation_optimization1=='Force (window size)'",
        numericInput('simulations_force_window',"Minimum (possible) window size to consider [bp]",min=100000,max=59128983,value=1000000)
    )})
    
    output$filtration_thresholdUI<-renderUI({conditionalPanel(
        condition="input.final_filter=='Yes'",
        numericInput('quality_filter','Quality Threshold',value=10.76,min=0,max=600,step = 0.01)
    )
    })
    
    #Analysis
    observeEvent(input$do,{
        shinyjs::html("text", paste0("<br>Starting analysis with CopyDetective...<br><br>"), add = FALSE)
        if(file.exists(input$input_folder)==F){
            shinyjs::html("text", paste0("Input folder does not exist","<br>"), add = TRUE) 
            return()
        }
        if(file.exists(input$output_folder)==F){
            shinyjs::html("text", paste0("Output folder does not exist","<br>"), add = TRUE) 
            return()
        }
        dir2<-paste0(input$input_folder,"/")
        dir_out<-paste0(input$output_folder,"/")
        
        samples_temp<-input$sampleFile
        if(is.null(samples_temp)){
            shinyjs::html("text", paste0("No sample names file defined","<br>"), add = TRUE) 
            return()
        }
        samples<-read.table(samples_temp$datapath,header=F,quote = "",sep="\t",stringsAsFactors = F)
        if(sum(is.na(samples[,1]))!=0||sum(is.na(samples[,2]))!=0){
            shinyjs::html("text", paste0("At least one sample without matching 2nd sample required","<br>"), add = TRUE) 
            return()
        }
        samples_g<-data.frame(Germline=samples[,1],stringsAsFactors=F)
        samples_t<-data.frame(Tumor=samples[,2],stringsAsFactors=F)

        genome<-c(249250621,
                  492449994,
                  690472424,
                  881626700,
                  1062541960,
                  1233657027,
                  1392795690,
                  1539159712,
                  1680373143,
                  1815907890,
                  1950914406,
                  2084766301,
                  2199936179,
                  2307285719,
                  2409817111,
                  2500171864,
                  2581367074,
                  2659444322,
                  2718573305,
                  2781598825,
                  2829728720,
                  2881033286)
        
        
        detectionThresholds<-data.frame(Sample=samples_t[,1],Window_del=NA,Rate_del=NA,
                                        Window_dup=NA,Rate_dup=NA)
        
        erwartungswerte<-data.frame(Sample=samples_t[,1],Germline=NA)

        for(n in 1:length(samples_t[,1])){
            shinyjs::html("text", paste0("<br>Starting analysis with CopyDetective...<br><br>"), add = FALSE)
            shinyjs::html("text", paste0("Analyzing sample ",samples_t[n,1],"<br>"), add = TRUE)
            progress_sample <- shiny::Progress$new()
            progress_sample$set(message = paste0("Analyzing sample ",samples_t[n,1]), value = 0)
            output$sample<-renderText({paste0("Sample ",samples_t[n,1])})
            output$sample2<-renderText({paste0("Sample ",samples_t[n,1])})
            
            germline<-read.table(paste0(input$input_folder,"/",samples_g[n,1],".txt"),header = T,sep="\t",stringsAsFactors = F)
            tumor<-read.table(paste0(input$input_folder,"/",samples_t[n,1],".txt"),header = T,sep="\t",stringsAsFactors = F)
            germline<-cbind(germline[,c(1:4)],Sample=samples_g[n,1],germline[,c(5:8)],stringsAsFactors = F)
            tumor<-cbind(tumor[,c(1:4)],Sample=samples_t[n,1],tumor[,c(5:8)],stringsAsFactors = F)
            
            tumor<-tumor[order(tumor[,1],tumor[,2]),]
            germline<-germline[order(germline[,1],germline[,2]),]
            
            erwartungswerte[n,2]<-mean(germline[,9])
            
            #Preparations: determine the cell fractions and the confidence intervals + information on coverage and the strand
            #for tumor
            progress_preprocess <- shiny::Progress$new()
            progress_preprocess$set(message = paste0("Prepare case sample"), value = 0)
            del<-data.frame(upper=NA,mean=NA,lower=NA,cov=NA,flag=0,cov2=NA)
            dup<-data.frame(upper=NA,mean=NA,lower=NA,cov=NA,flag=0,cov2=NA)
            for(i in 1:length(tumor[,1])){
                progress_preprocess$inc(1/length(tumor[,1]))
                if(!is.na(tumor[i,9])){
                    lower<-binom.test(x=tumor[i,7],tumor[i,8],tumor[i,9],conf.level = 0.95)$conf.int[1]
                    upper<-binom.test(x=tumor[i,7],tumor[i,8],tumor[i,9],conf.level = 0.95)$conf.int[2]
                    if(tumor[i,9]<=0.5){
                        del[i,1]<-(2*lower-1)/(lower-1)
                        del[i,2]<-(2*tumor[i,9]-1)/(tumor[i,9]-1)
                        del[i,3]<-(2*upper-1)/(upper-1)
                        del[i,4]<-tumor[i,8]<mean(tumor[,8],na.rm=T)
                        del[i,5]<-0
                        del[i,6]<-tumor[i,8]
                        
                        if(upper>=1/3&&tumor[i,9]!=0){
                            dup[i,1]<-(2*lower-1)/(-1*lower)
                            dup[i,2]<-(2*tumor[i,9]-1)/(-1*tumor[i,9])
                            dup[i,3]<-(2*upper-1)/(-1*upper)
                            dup[i,4]<-tumor[i,8]<mean(tumor[,8],na.rm=T)
                            dup[i,5]<-0
                            dup[i,6]<-tumor[i,8]
                        }
                    }
                    
                    if(tumor[i,9]>0.5){
                        del[i,3]<-(2*lower-1)/(lower)
                        del[i,2]<-(2*tumor[i,9]-1)/(tumor[i,9])
                        del[i,1]<-(2*upper-1)/(upper)
                        del[i,4]<-tumor[i,8]<mean(tumor[,8],na.rm=T)
                        del[i,5]<-1
                        del[i,6]<-tumor[i,8]
                        
                        if(lower<=2/3&&tumor[i,9]!=1){
                            dup[i,3]<-(2*lower-1)/(1-lower)
                            dup[i,2]<-(2*tumor[i,9]-1)/(1-tumor[i,9])
                            dup[i,1]<-(2*upper-1)/(1-upper)
                            dup[i,4]<-tumor[i,8]<mean(tumor[,8],na.rm=T)
                            dup[i,5]<-1
                            dup[i,6]<-tumor[i,8]
                        }
                    }
                }
            }
            tumor[is.na(tumor[,8]),8]<-1
            progress_preprocess$close()
            
            #for germline
            progress_preprocess <- shiny::Progress$new()
            progress_preprocess$set(message = paste0("Prepare control sample"), value = 0)
            del_g<-data.frame(upper=NA,mean=NA,lower=NA,cov=NA,flag=0,cov2=NA)
            dup_g<-data.frame(upper=NA,mean=NA,lower=NA,cov=NA,flag=0,cov2=NA)
            for(i in 1:length(tumor[,1])){
                #progress_preprocess$inc(1/length(tumor[,1]))
                if(!is.na(germline[i,9])){
                    lower<-binom.test(x=germline[i,7],germline[i,8],germline[i,9],conf.level = 0.95)$conf.int[1]
                    upper<-binom.test(x=germline[i,7],germline[i,8],germline[i,9],conf.level = 0.95)$conf.int[2]
                    if(germline[i,9]<=0.5){
                        del_g[i,1]<-(2*lower-1)/(lower-1)
                        del_g[i,2]<-(2*germline[i,9]-1)/(germline[i,9]-1)
                        del_g[i,3]<-(2*upper-1)/(upper-1)
                        del_g[i,4]<-germline[i,8]<mean(germline[,8],na.rm=T)
                        del_g[i,5]<-0
                        del_g[i,6]<-germline[i,8]
                        
                        if(upper>=1/3&&germline[i,9]!=0){
                            dup_g[i,1]<-(2*lower-1)/(-1*lower)
                            dup_g[i,2]<-(2*germline[i,9]-1)/(-1*germline[i,9])
                            dup_g[i,3]<-(2*upper-1)/(-1*upper)
                            dup_g[i,4]<-germline[i,8]<mean(germline[,8],na.rm=T)
                            dup_g[i,5]<-0
                            dup_g[i,6]<-germline[i,8]
                        }
                    }
                    
                    if(germline[i,9]>0.5){
                        del_g[i,3]<-(2*lower-1)/(lower)
                        del_g[i,2]<-(2*germline[i,9]-1)/(germline[i,9])
                        del_g[i,1]<-(2*upper-1)/(upper)
                        del_g[i,4]<-germline[i,8]<mean(germline[,8],na.rm=T)
                        del_g[i,5]<-1
                        del_g[i,6]<-germline[i,8]
                        
                        if(lower<=2/3&&germline[i,9]!=1){
                            dup_g[i,3]<-(2*lower-1)/(1-lower)
                            dup_g[i,2]<-(2*germline[i,9]-1)/(1-germline[i,9])
                            dup_g[i,1]<-(2*upper-1)/(1-upper)
                            dup_g[i,4]<-germline[i,8]<mean(germline[,8],na.rm=T)
                            dup_g[i,5]<-1
                            dup_g[i,6]<-germline[i,8]
                        }
                    }
                    if(germline[i,9]==0.5){
                        del_g[i,2]<-dup_g[i,2]<-0
                    }
                }
            }
            progress_preprocess$close()
            
            progress_sample$inc(1/4)
            if(input$simulation=="No"){
                if(input$exactSwitch=="Simulation"){
                    shinyjs::html("text", paste0("<br>","&nbsp&nbsp&nbspDetermining thresholds for CNV detection (simulation approach)","<br>"), add = TRUE)
                    rate<-seq(input$simulate_rates[2],input$simulate_rates[1],(-1)*input$simulate_steps)
                    rate<-rate/100
                    thresholds_del<-data.frame(rate=rate,minSnps=NA,window=NA)
                    thresholds_dup<-data.frame(rate=rate,minSnps=NA,window=NA)
                    shinyjs::html("text", paste0("&nbsp&nbsp&nbsp&nbsp&nbsp&nbspA: Evaluating Cell Fractions","<br>"), add = TRUE)  
                    progress <- shiny::Progress$new()
                    progress$set(message = "A: Evaluating Cell Fractions", value = 0)
                    remember_window<-0
                    for(j in 1:length(rate)){
                        #message("Rate: ",rate[j])
                        progress$inc(input$simulate_steps/100)
                        thresholds_snps<-data.frame(snps=1:100,significant_del=0,significant_dup=0)
                        vaf_del<-(rate[j]-1)/(-2+rate[j])
                        vaf_g<-0.5  
                        vaf_dup<-1/(2+rate[j])
                        simulations<-input$simulations_nr
                        
                        info_del<-data.frame(cov_tumor=rep(NA,100*simulations),vaf_tumor_del=rep(NA,100*simulations),
                                             vaf_tumor_lower_del=rep(NA,100*simulations),vaf_tumor_upper_del=rep(NA,100*simulations),
                                             cov_germline=rep(NA,100*simulations),vaf_germline=rep(NA,100*simulations),
                                             vaf_germline_lower=rep(NA,100*simulations),vaf_germline_upper=rep(NA,100*simulations),
                                             vaf_tumor_dup=rep(NA,100*simulations),
                                             vaf_tumor_lower_dup=rep(NA,100*simulations),vaf_tumor_upper_dup=rep(NA,100*simulations))
                        tumor[which(tumor[,8]==0),8]<-NA
                        temp<-round(rlnorm(meanlog=mean(log(tumor[,8]),na.rm=T),sdlog = sd(log(tumor[,8]),na.rm=T),n=200*simulations))
                        info_del[,1]<-temp[temp>0][1:(100*simulations)]
                        temp<-round(rlnorm(meanlog=mean(log(germline[,8])),sdlog = sd(log(germline[,8])),n=200*simulations))
                        info_del[,5]<-temp[temp>0][1:(100*simulations)]
                        
                        info_del[,2]<-apply(info_del,MARGIN=1,
                                            FUN=function(x){sqrt((round(rnorm(mean=vaf_del*x[1],
                                                                              sd = 0.05*vaf_del*x[1],n=1)))**2)})/info_del[,1]
                        info_del[,6]<-apply(info_del,MARGIN=1,
                                            FUN=function(x){sqrt((round(rnorm(mean=vaf_g*x[5],
                                                                              sd = 0.05*vaf_g*x[5],n=1)))**2)})/info_del[,5]
                        info_del[,9]<-apply(info_del,MARGIN=1,
                                            FUN=function(x){sqrt((round(rnorm(mean=vaf_dup*x[1],
                                                                              sd = 0.05*vaf_dup*x[1],n=1)))**2)})/info_del[,1]
                        
                        info_del_r<-data.frame(cov_tumor=info_del[,1],r_tumor_del=rep(NA,100*simulations),
                                               r_tumor_lower_del=rep(NA,100*simulations),r_tumor_upper_del=rep(NA,100*simulations),
                                               cov_germline=info_del[,5],r_germline_del=rep(NA,100*simulations),
                                               r_germline_del_lower=rep(NA,100*simulations),r_germline_del_upper=rep(NA,100*simulations),
                                               r_tumor_dup=rep(NA,100*simulations),
                                               r_tumor_lower_dup=rep(NA,100*simulations),r_tumor_upper_dup=rep(NA,100*simulations),
                                               r_germline_dup=rep(NA,100*simulations),
                                               r_germline_dup_lower=rep(NA,100*simulations),r_germline_dup_upper=rep(NA,100*simulations))
                        
                        #deletions
                        #tumor
                        for(k in c(2)){
                            temp1<-info_del[,k]
                            temp1[info_del[,2]>0.5]<-NA
                            temp1.1<-(2*temp1-1)/(temp1-1)
                            temp2<-info_del[,k]
                            temp2[info_del[,2]<=0.5]<-NA
                            temp2.1<-(2*temp2-1)/temp2
                            temp3<-apply(cbind(temp1.1,temp2.1),MARGIN=1,FUN=function(x){max(x,na.rm=T)})
                            info_del_r[,k]<-temp3
                        }
                        #germline
                        for(k in c(6)){
                            temp1<-info_del[,k]
                            temp1[info_del[,6]>0.5]<-NA
                            temp1.1<-(2*temp1-1)/(temp1-1)
                            temp2<-info_del[,k]
                            temp2[info_del[,6]<=0.5]<-NA
                            temp2.1<-(2*temp2-1)/temp2
                            temp3<-apply(cbind(temp1.1,temp2.1),MARGIN=1,FUN=function(x){max(x,na.rm=T)})
                            info_del_r[,k]<-temp3
                        }
                        #duplications
                        #tumor
                        for(k in c(9)){
                            temp1<-info_del[,k]
                            temp1[info_del[,9]>0.5]<-NA
                            temp1.1<-(2*temp1-1)/(-1*temp1)
                            temp2<-info_del[,k]
                            temp2[info_del[,9]<=0.5]<-NA
                            temp2.1<-(2*temp2-1)/(1-temp2)
                            temp3<-apply(cbind(temp1.1,temp2.1),MARGIN=1,FUN=function(x){max(x,na.rm=T)})
                            info_del_r[,k]<-temp3
                        }
                        #germline
                        for(k in c(6)){
                            temp1<-info_del[,k]
                            temp1[info_del[,6]>0.5]<-NA
                            temp1.1<-(2*temp1-1)/(-1*temp1)
                            temp2<-info_del[,k]
                            temp2[info_del[,6]<=0.5]<-NA
                            temp2.1<-(2*temp2-1)/(1-temp2)
                            temp3<-apply(cbind(temp1.1,temp2.1),MARGIN=1,FUN=function(x){max(x,na.rm=T)})
                            info_del_r[,k+6]<-temp3
                        }
                        
                        test_del2<-cbind(info_del_r[,2],info_del_r[,6],info_del_r[,1],info_del_r[,5])
                        test_dup2<-cbind(info_del_r[,9],info_del_r[,12],info_del_r[,1],info_del_r[,5])
                        
                        if(remember_window==0){
                            i<-2   
                        }
                        if(remember_window!=0){
                            i<-remember_window
                        }
                        del_enough<-F
                        dup_enough<-F
                        start<-seq(1,simulations*100,100)
                        while(i<=100&&(del_enough==F||dup_enough==F)){
                            for(k in 1:length(start)){
                                if(del_enough==F){
                                    temp<-wtd.t.test(x=test_del2[start[k]:(start[k]+i-1),1],
                                                     y=test_del2[start[k]:(start[k]+i-1),2],
                                                     weight=test_del2[start[k]:(start[k]+i-1),3],
                                                     weighty=test_del2[start[k]:(start[k]+i-1),4],
                                                     alternative="greater")
                                    if(!is.na(temp$coefficients[3])&&temp$coefficients[3]<0.0125){
                                        thresholds_snps[i,2]<-sum(thresholds_snps[i,2],1)
                                    } 
                                }
                                if(dup_enough==F){
                                    temp<-wtd.t.test(x=test_dup2[start[k]:(start[k]+i-1),1],
                                                     y=test_dup2[start[k]:(start[k]+i-1),2],
                                                     weight=test_dup2[start[k]:(start[k]+i-1),3],
                                                     weighty=test_dup2[start[k]:(start[k]+i-1),4],
                                                     alternative="greater")
                                    if(!is.na(temp$coefficients[3])&&temp$coefficients[3]<0.0125){
                                        thresholds_snps[i,3]<-sum(thresholds_snps[i,3],1)
                                    }   
                                }
                            }
                            if(thresholds_snps[i,2]>=simulations){
                                del_enough<-T
                            }
                            if(thresholds_snps[i,3]>=simulations){
                                dup_enough<-T
                            }
                            i<-i+1
                        }
                        thresholds_del[rate[j]==thresholds_del[,1],2]<-thresholds_snps[min(which(thresholds_snps[,2]>=input$aimSens*simulations),100,na.rm=T),1]
                        thresholds_dup[rate[j]==thresholds_dup[,1],2]<-thresholds_snps[min(which(thresholds_snps[,3]>=input$aimSens*simulations),100,na.rm=T),1]
                        remember_window<-min(thresholds_del[rate[j]==thresholds_del[,1],2],thresholds_dup[rate[j]==thresholds_dup[,1],2])
                    }
                    
                    thresholds_del<-thresholds_del[thresholds_del[,2]<100,]
                    thresholds_dup<-thresholds_dup[thresholds_dup[,2]<100,]
                    progress$close()
                }
                if(input$exactSwitch=="Exact"){
                    shinyjs::html("text", paste0("<br>","&nbsp&nbsp&nbspDetermining thresholds for CNV detection (exact approach)","<br>"), add = TRUE)
                    rate<-seq(input$simulate_rates[2],input$simulate_rates[1],(-1)*input$simulate_steps)
                    rate<-rate/100
                    thresholds_del<-data.frame(rate=rate,minSnps=NA,window=NA)
                    thresholds_dup<-data.frame(rate=rate,minSnps=NA,window=NA)
                    shinyjs::html("text", paste0("&nbsp&nbsp&nbsp&nbsp&nbsp&nbspA: Evaluating Cell Fractions","<br>"), add = TRUE)  
                    progress <- shiny::Progress$new()
                    progress$set(message = "A: Evaluating Cell Fractions", value = 0)
                    remember_window<-0
                    for(j in 1:length(rate)){
                        #message("Rate: ",rate[j])
                        progress$inc(input$simulate_steps/100)
                        thresholds_snps<-data.frame(snps=1:100,significant_del=0,significant_dup=0)
                        vaf_del<-(rate[j]-1)/(-2+rate[j])
                        vaf_g<-0.5  
                        vaf_dup<-1/(2+rate[j])

                        info_del<-data.frame(cov_tumor=tumor[,8],vaf_tumor_del=NA,
                                             vaf_tumor_lower_del=NA,vaf_tumor_upper_del=NA,
                                             cov_germline=germline[,8],vaf_germline=germline[,9],
                                             vaf_germline_lower=NA,vaf_germline_upper=NA,
                                             vaf_tumor_dup=NA,
                                             vaf_tumor_lower_dup=NA,NA)
                        info_del[which(info_del[,1]==0),1]<-NA

                        info_del[,2]<-apply(info_del,MARGIN=1,
                                            FUN=function(x){sqrt((round(rnorm(mean=vaf_del*x[1],
                                                                              sd = 0.05*vaf_del*x[1],n=1)))**2)})/info_del[,1]
                        info_del[,9]<-apply(info_del,MARGIN=1,
                                            FUN=function(x){sqrt((round(rnorm(mean=vaf_dup*x[1],
                                                                              sd = 0.05*vaf_dup*x[1],n=1)))**2)})/info_del[,1]
                        
                        if(input$simulate_germline_vaf=="Simulated (expected value 0.5)"){
                          info_del[,6]<-apply(info_del,MARGIN=1,
                                             FUN=function(x){sqrt((round(rnorm(mean=vaf_g*x[5],
                                                                               sd = 0.05*vaf_g*x[5],n=1)))**2)})/info_del[,5]
                        }
                        
                        
                        info_del_r<-data.frame(cov_tumor=info_del[,1],r_tumor_del=NA,
                                               r_tumor_lower_del=NA,r_tumor_upper_del=NA,
                                               cov_germline=info_del[,5],r_germline_del=NA,
                                               r_germline_del_lower=NA,r_germline_del_upper=NA,
                                               r_tumor_dup=NA,
                                               r_tumor_lower_dup=NA,r_tumor_upper_dup=NA,
                                               r_germline_dup=NA,
                                               r_germline_dup_lower=NA,r_germline_dup_upper=NA)
                        
                        #deletions
                        #tumor
                        for(k in c(2)){
                            temp1<-info_del[,k]
                            temp1[info_del[,2]>0.5]<-NA
                            temp1.1<-(2*temp1-1)/(temp1-1)
                            temp2<-info_del[,k]
                            temp2[info_del[,2]<=0.5]<-NA
                            temp2.1<-(2*temp2-1)/temp2
                            temp3<-apply(cbind(temp1.1,temp2.1),MARGIN=1,FUN=function(x){max(x,na.rm=T)})
                            info_del_r[,k]<-temp3
                        }
                        #germline
                        for(k in c(6)){
                            temp1<-info_del[,k]
                            temp1[info_del[,6]>0.5]<-NA
                            temp1.1<-(2*temp1-1)/(temp1-1)
                            temp2<-info_del[,k]
                            temp2[info_del[,6]<=0.5]<-NA
                            temp2.1<-(2*temp2-1)/temp2
                            temp3<-apply(cbind(temp1.1,temp2.1),MARGIN=1,FUN=function(x){max(x,na.rm=T)})
                            info_del_r[,k]<-temp3
                        }
                        #duplications
                        #tumor
                        for(k in c(9)){
                            temp1<-info_del[,k]
                            temp1[info_del[,9]>0.5]<-NA
                            temp1.1<-(2*temp1-1)/(-1*temp1)
                            temp2<-info_del[,k]
                            temp2[info_del[,9]<=0.5]<-NA
                            temp2.1<-(2*temp2-1)/(1-temp2)
                            temp3<-apply(cbind(temp1.1,temp2.1),MARGIN=1,FUN=function(x){max(x,na.rm=T)})
                            info_del_r[,k]<-temp3
                        }
                        #germline
                        for(k in c(6)){
                            temp1<-info_del[,k]
                            temp1[info_del[,6]>0.5]<-NA
                            temp1.1<-(2*temp1-1)/(-1*temp1)
                            temp2<-info_del[,k]
                            temp2[info_del[,6]<=0.5]<-NA
                            temp2.1<-(2*temp2-1)/(1-temp2)
                            temp3<-apply(cbind(temp1.1,temp2.1),MARGIN=1,FUN=function(x){max(x,na.rm=T)})
                            info_del_r[,k+6]<-temp3
                        }
                        
                        test_del2<-cbind(info_del_r[,2],info_del_r[,6],info_del_r[,1],info_del_r[,5])
                        test_dup2<-cbind(info_del_r[,9],info_del_r[,12],info_del_r[,1],info_del_r[,5])
                        
                        if(remember_window==0){
                            test_window<-2   
                        }
                        if(remember_window!=0){
                            test_window<-remember_window
                        }
                        del_enough<-F
                        dup_enough<-F
                        while(test_window<=100&&(del_enough==F||dup_enough==F)){
                            #message("Rate: ",j," window: ",test_window)
                            for(k in 1:(length(test_del2[,1])-test_window+1)){
                                if(del_enough==F){
                                    temp<-wtd.t.test(x=test_del2[k:(k+test_window-1),1],
                                                     y=test_del2[k:(k+test_window-1),2],
                                                     weight=test_del2[k:(k+test_window-1),3],
                                                     weighty=test_del2[k:(k+test_window-1),4],
                                                     alternative="greater")
                                    if(!is.na(temp$coefficients[3])&&temp$coefficients[3]<0.0125){
                                        thresholds_snps[test_window,2]<-sum(thresholds_snps[test_window,2],1)
                                    } 
                                }
                                if(dup_enough==F){
                                    temp<-wtd.t.test(x=test_dup2[k:(k+test_window-1),1],
                                                     y=test_dup2[k:(k+test_window-1),2],
                                                     weight=test_dup2[k:(k+test_window-1),3],
                                                     weighty=test_dup2[k:(k+test_window-1),4],
                                                     alternative="greater")
                                    if(!is.na(temp$coefficients[3])&&temp$coefficients[3]<0.0125){
                                        thresholds_snps[test_window,3]<-sum(thresholds_snps[test_window,3],1)
                                    }   
                                }
                            }
                            if(thresholds_snps[test_window,2]>=input$aimSens*(length(test_del2[,1])-test_window+1)){
                                del_enough<-T
                            }
                            if(thresholds_snps[test_window,3]>=input$aimSens*(length(test_del2[,1])-test_window+1)){
                                dup_enough<-T
                            }
                            test_window<-test_window+1
                        }
                        thresholds_del[rate[j]==thresholds_del[,1],2]<-thresholds_snps[min(which(thresholds_snps[,2]>=input$aimSens*(length(test_del2[,1])-test_window+1)),100,na.rm=T),1]
                        thresholds_dup[rate[j]==thresholds_dup[,1],2]<-thresholds_snps[min(which(thresholds_snps[,3]>=input$aimSens*(length(test_del2[,1])-test_window+1)),100,na.rm=T),1]
                        remember_window<-min(thresholds_del[rate[j]==thresholds_del[,1],2],thresholds_dup[rate[j]==thresholds_dup[,1],2])
                    }
                    
                    thresholds_del<-thresholds_del[thresholds_del[,2]<100,]
                    thresholds_dup<-thresholds_dup[thresholds_dup[,2]<100,]
                    progress$close()
                }
                

                
                shinyjs::html("text", paste0("&nbsp&nbsp&nbsp&nbsp&nbsp&nbspB: Evaluating Window Sizes","<br>"), add = TRUE)  
                progress <- shiny::Progress$new()
                progress$set(message = "B: Evaluating Window Sizes for deletions", value = 0)
                snps_to_consider<-unique(thresholds_del[,2])
                for(j in snps_to_consider){
                    progress$inc(1/(length(snps_to_consider)+1))
                    window<-data.frame(chr=NA,start=NA,end=NA)
                    counter<-1
                    for(chr in 1:22){
                        tumor_1<-tumor[tumor[,1]==chr,]
                        if((length(tumor_1[,1])-j)>=1){
                            for(i in 1:(length(tumor_1[,1])-j)){
                                window[counter,1]<-tumor_1[i,1]
                                window[counter,2]<-tumor_1[i,2]
                                window[counter,3]<-tumor_1[(i+j),2]
                                counter<-counter+1
                            }
                        }
                        if((length(tumor_1[,1])-j)<1){
                            window[counter,1]<-tumor_1[1,1]
                            window[counter,2]<-1
                            if(tumor[1,1]=="1"){
                                window[counter,3]<-genome[1]
                            }
                            if(tumor[1,1]!="1"){
                                window[counter,3]<-genome[as.numeric(tumor_1[1,1])]-genome[as.numeric(tumor_1[1,1])-1]
                            }
                            counter<-counter+1
                        }
                    }  
                    window<-cbind(window,Length=(window[,3]-window[,2]))
                    thresholds_del[thresholds_del[,2]==j,3]<-quantile(window[,4],(input$window_percentile/100))
                }
                thresholds_del<-cbind(thresholds_del,windowTransformed=thresholds_del[,3]/max(thresholds_del[,3],na.rm=T))
                thresholds_del<-cbind(thresholds_del,Distance=sqrt(thresholds_del[,1]**2+thresholds_del[,4]**2))
                progress$close()
                
                progress <- shiny::Progress$new()
                progress$set(message = "B: Evaluating Window Sizes for duplications", value = 0)
                snps_to_consider<-unique(thresholds_dup[,2])
                for(j in snps_to_consider){
                    progress$inc(1/(length(snps_to_consider)+1))
                    tumor_temp<-tumor[!is.na(dup[,1]),]
                    window<-data.frame(chr=NA,start=NA,end=NA)
                    counter<-1
                    for(chr in 1:22){
                        tumor_1<-tumor_temp[tumor_temp[,1]==chr,]
                        if((length(tumor_1[,1])-j)>=1){
                            for(i in 1:(length(tumor_1[,1])-j)){
                                window[counter,1]<-tumor_1[i,1]
                                window[counter,2]<-tumor_1[i,2]
                                window[counter,3]<-tumor_1[(i+j),2]
                                counter<-counter+1
                            }
                        }
                        if((length(tumor_1[,1])-j)<1){
                            window[counter,1]<-tumor_1[1,1]
                            window[counter,2]<-1
                            if(tumor[1,1]=="1"){
                                window[counter,3]<-genome[1]
                            }
                            if(tumor[1,1]!="1"){
                                window[counter,3]<-genome[as.numeric(tumor_1[1,1])]-genome[as.numeric(tumor_1[1,1])-1]
                            }
                            counter<-counter+1
                        }
                    }  
                    window<-cbind(window,Length=(window[,3]-window[,2]))
                    thresholds_dup[thresholds_dup[,2]==j,3]<-quantile(window[,4],(input$window_percentile/100))
                }
                thresholds_dup<-cbind(thresholds_dup,windowTransformed=thresholds_dup[,3]/max(thresholds_dup[,3],na.rm=T))
                thresholds_dup<-cbind(thresholds_dup,Distance=sqrt(thresholds_dup[,1]**2+thresholds_dup[,4]**2)) 
                
                progress$close()
                
                
                shinyjs::html("text", paste0("<br>","&nbsp&nbsp&nbspDetection thresholds for CNVs:","<br>"), add = TRUE)
                if((!is.na(input$simulation_optimization)&&input$simulation_optimization=="Compromize (cell fraction and window size)")||
                   (!is.na(input$simulation_optimization1)&&input$simulation_optimization1=="Compromize (cell fraction and window size)")){
                    del_rate<-thresholds_del[which.min(thresholds_del[,5]),1]
                    del_window<-thresholds_del[which.min(thresholds_del[,5]),3]
                    del_snps<-thresholds_del[which.min(thresholds_del[,5]),2]
                    shinyjs::html("text", paste0("&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbspDeletions: ",del_window,
                                                 " bp; ",del_rate*100,"% cells","<br>"), add = TRUE)
                    
                    dup_rate<-thresholds_dup[which.min(thresholds_dup[,5]),1]
                    dup_window<-thresholds_dup[which.min(thresholds_dup[,5]),3]
                    dup_snps<-thresholds_dup[which.min(thresholds_dup[,5]),2]
                    shinyjs::html("text", paste0("&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbspDuplications: ",dup_window,
                                                 " bp; ",dup_rate*100,"% cells","<br>"), add = TRUE)
                }
                if((!is.na(input$simulation_optimization)&&input$simulation_optimization=="Force (cell fraction)")||
                   (!is.na(input$simulation_optimization1)&&input$simulation_optimization1=="Force (cell fraction)")){
                    if(input$simulations_force_rate<=100*min(thresholds_del[!is.na(thresholds_del[,5]),1])){
                        goforit<-which.min(thresholds_del[!is.na(thresholds_del[,5]),1])
                    }
                    if(input$simulations_force_rate>100*min(thresholds_del[!is.na(thresholds_del[,5]),1])){
                        goforit<-which(thresholds_del[!is.na(thresholds_del[,5]),1]==(input$simulations_force_rate)/100)                       
                    }
                    del_rate<-thresholds_del[goforit,1]
                    del_window<-thresholds_del[goforit,3]
                    del_snps<-thresholds_del[goforit,2]
                    shinyjs::html("text", paste0("&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbspDeletions: ",del_window,
                                                 " bp; ",del_rate*100,"% cells","<br>"), add = TRUE) 
                    
                    if(input$simulations_force_rate<=100*min(thresholds_dup[!is.na(thresholds_dup[,5]),1])){
                        goforit<-which.min(thresholds_dup[!is.na(thresholds_dup[,5]),1])
                    }
                    if(input$simulations_force_rate>100*min(thresholds_dup[!is.na(thresholds_dup[,5]),1])){
                        goforit<-which(thresholds_dup[!is.na(thresholds_dup[,5]),1]==(input$simulations_force_rate)/100)                       
                    }
                    dup_rate<-thresholds_dup[goforit,1]
                    dup_window<-thresholds_dup[goforit,3]
                    dup_snps<-thresholds_dup[goforit,2]
                    shinyjs::html("text", paste0("&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbspDuplications: ",dup_window,
                                                 " bp; ",dup_rate*100,"% cells","<br>"), add = TRUE)
                }
                if((!is.na(input$simulation_optimization)&&input$simulation_optimization=="Force (window size)")||
                   (!is.na(input$simulation_optimization1)&&input$simulation_optimization1=="Force (window size)")){
                    if(input$simulations_force_window<=min(thresholds_del[!is.na(thresholds_del[,5]),3])){
                        goforit<-which.min(thresholds_del[!is.na(thresholds_del[,5]),3])
                    }
                    if(input$simulations_force_window>min(thresholds_del[!is.na(thresholds_del[,5]),3])){
                        goforit<-which(thresholds_del[!is.na(thresholds_del[,5]),3]==input$simulations_force_window)                       
                    }
                    del_rate<-thresholds_del[goforit,1]
                    del_window<-thresholds_del[goforit,3]
                    del_snps<-thresholds_del[goforit,2]
                    shinyjs::html("text", paste0("&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbspDeletions: ",del_window,
                                                 " bp; ",del_rate*100,"% cells","<br>"), add = TRUE) 
                    
                    if(input$simulations_force_window<=min(thresholds_dup[!is.na(thresholds_dup[,5]),3])){
                        goforit<-which.min(thresholds_dup[!is.na(thresholds_dup[,5]),3])
                    }
                    if(input$simulations_force_window>min(thresholds_dup[!is.na(thresholds_dup[,5]),3])){
                        goforit<-which(thresholds_dup[!is.na(thresholds_dup[,5]),3]==input$simulations_force_window)                       
                    }
                    dup_rate<-thresholds_dup[goforit,1]
                    dup_window<-thresholds_dup[goforit,3]
                    dup_snps<-thresholds_dup[goforit,2]
                    shinyjs::html("text", paste0("&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbspDuplications: ",dup_window,
                                                 " bp; ",dup_rate*100,"% cells","<br>"), add = TRUE)
                }
                
                if(n==1){
                    detectionThresholds<-data.frame(Sample=c(samples_t[n,1]),Window_del=c(del_window),
                                                    CF_del=c(del_rate),SNPs_del=c(del_snps),Window_dup=c(dup_window),
                                                    CF_dup=c(dup_rate),SNPs_dup=c(dup_snps))
                }
                if(n>1){
                    detectionThresholds<-rbind(detectionThresholds,data.frame(Sample=c(samples_t[n,1]),
                                                                              Window_del=c(del_window),
                                                                              CF_del=c(del_rate),
                                                                              SNPs_del=c(del_snps),
                                                                              Window_dup=c(dup_window),
                                                                              CF_dup=c(dup_rate),
                                                                              SNPs_dup=c(dup_snps)))
                }
                output$table_dt <- renderDataTable(datatable(detectionThresholds))
                write.table(detectionThresholds,paste0(input$output_folder,"/DetectionThresholds.txt"),
                            sep="\t",row.names = F,quote=F)
            }
            if(input$simulation=="Yes"){
                detection_temp<-input$detectionThresholdsFile
                detectionThresholds<-read.table(detection_temp$datapath,header=T,quote = "",sep="\t",stringsAsFactors = F)
                output$table_dt <- renderDataTable(datatable(detectionThresholds))
            }

            shinyjs::html("text", paste0("<br>","&nbsp&nbsp&nbspDetecting CNVs","<br>"), add = TRUE)
            progress_sample$inc(1/4)
            
            #Calling deletions:
            window<-detectionThresholds[n,2]
            results<-data.frame(chr=NA,start=NA,end=NA,variant=NA,p.value=NA,cells=NA,
                                snps=NA,qual=NA,sd=NA)
            keep_in<-rep(F,length(del[,1]))
            counter<-1
            progress <- shiny::Progress$new()
            progress$set(message = "Detecting deletions", value = 0)
            for(chr in 1:22){
                progress$inc(1/23)
                tumor_1<-tumor[tumor[,1]==chr,]
                del_1<-del[tumor[,1]==chr,]
                del_g1<-del_g[tumor[,1]==chr,]
                which_chr<-which(tumor[,1]==chr)
                for(i in 0:(length(tumor_1[,1])-1)){
                    if(i==0){
                        start<-1
                    }
                    if(i>0){
                        start<-tumor_1[i,2]+1
                    }
                    end<-start+window
                    
                    del_oi<-del_1[tumor_1[,2]>=start&tumor_1[,2]<=(start+window),]
                    del_goi<-del_g1[tumor_1[,2]>=start&tumor_1[,2]<=(start+window),]
                    if(length(del_oi[,1])>0&&length(del_goi[,1])>0&&
                       sum(!is.na(del_oi[,2]))>1&&sum(!is.na(del_goi[,2]))>1&&
                       sum(!is.na(del_oi[,2]-del_goi[,2]))>1&&
                       sum(!is.na(unique(round(del_oi[,2]-del_goi[,2],digits = 7))))>1){
                        not_na<-intersect(which(!is.na(del_oi[,1])),which(!is.na(del_goi[,1])))
                        del_oi<-del_oi[not_na,]
                        del_goi<-del_goi[not_na,]
                        temp<-wtd.t.test(x=del_oi[,2],y=del_goi[,2],weight=del_oi[,6],weighty=del_goi[,6],alternative = "greater")
                        temp2<-wtd.t.test(x=del_oi[,2],y=del_goi[,2],weight=del_oi[,6],weighty=del_goi[,6],alternative = "less")

                        results[counter,1]<-chr
                        results[counter,2]<-start
                        results[counter,3]<-start+window
                        results[counter,4]<-"del"
                        results[counter,7]<-length(del_goi[,1])
                        #CNV in more tumor cells
                        if(temp$additional[1]>0){
                            results[counter,5]<-temp$coefficients[3]
                            results[counter,6]<-temp$additional[1]
                            results[counter,8]<-temp$coefficients[3]**(-1)
                            results[counter,9]<-temp$additional[4]
                        }
                        #CNV in more germline cells -> use later
                        if(temp$additional[1]<0){
                            results[counter,5]<-temp2$coefficients[3]
                            results[counter,6]<-temp2$additional[1]
                            results[counter,8]<--1*temp2$coefficients[3]**(-1)
                            results[counter,9]<-temp2$additional[4]
                        }
                        if(temp$additional[1]==0){
                            results[counter,5]<-temp$coefficients[3]
                            results[counter,6]<-temp$additional[1]
                            results[counter,8]<-0
                            results[counter,9]<-temp$additional[4]
                        }
                        if(results[counter,8]>=80){
                            which_pos<-which(tumor[,2]>=start&tumor[,2]<=end)
                            which_inter<-intersect(which_chr,which_pos)
                            keep_in[which_inter]<-T
                        }
                        counter<-counter+1
                    }
                }
            } 
            del[keep_in==F,2]<-NA
            progress$close()
            
            #Calling duplications:
            keep_in<-rep(F,length(dup[,1]))
            window<-detectionThresholds[n,5]
            progress <- shiny::Progress$new()
            progress$set(message = "Detecting duplications", value = 0)
            for(chr in 1:22){
                progress$inc(1/23)
                tumor_1<-tumor[tumor[,1]==chr,]
                dup_1<-dup[tumor[,1]==chr,]
                dup_g1<-dup_g[tumor[,1]==chr,]
                which_chr<-which(tumor[,1]==chr)
                for(i in 0:(length(tumor_1[,1])-1)){
                    if(i==0){
                        start<-1
                    }
                    if(i>0){
                        start<-tumor_1[i,2]+1
                    }
                    end<-start+window
                    
                    dup_oi<-dup_1[tumor_1[,2]>=start&tumor_1[,2]<=(start+window),]
                    dup_goi<-dup_g1[tumor_1[,2]>=start&tumor_1[,2]<=(start+window),]
                    
                    if(length(dup_oi[,1])>0&&length(dup_goi[,1])>0&&
                       sum(!is.na(dup_oi[,2]))>1&&sum(!is.na(dup_goi[,2]))>1&&
                       sum(!is.na(dup_oi[,2]-dup_goi[,2]))>1&&
                       sum(!is.na(unique(round(dup_oi[,2]-dup_goi[,2],digits = 7))))>1){
                        not_na<-intersect(which(!is.na(dup_oi[,1])),which(!is.na(dup_goi[,1])))
                        dup_oi<-dup_oi[not_na,]
                        dup_goi<-dup_goi[not_na,]
                        temp<-wtd.t.test(x=dup_oi[,2],y=dup_goi[,2],weight=dup_oi[,6],weighty=dup_goi[,6],alternative = "greater")
                        temp2<-wtd.t.test(x=dup_oi[,2],y=dup_goi[,2],weight=dup_oi[,6],weighty=dup_goi[,6],alternative = "less")
                        results[counter,1]<-chr
                        results[counter,2]<-start
                        results[counter,3]<-start+window
                        results[counter,4]<-"dup"
                        results[counter,7]<-length(dup_goi[,1])
                        #CNV in more tumor cells
                        if(temp$additional[1]>0){
                            results[counter,5]<-temp$coefficients[3]
                            results[counter,6]<-temp$additional[1]
                            results[counter,8]<-temp$coefficients[3]**(-1)
                            results[counter,9]<-temp$additional[4]
                        }
                        #CNV in more germline cells -> use later
                        if(temp$additional[1]<0){
                            results[counter,5]<-temp2$coefficients[3]
                            results[counter,6]<-temp2$additional[1]
                            results[counter,8]<--1*temp2$coefficients[3]**(-1)
                            results[counter,9]<-temp2$additional[4]
                        }
                        if(temp$additional[1]==0){
                            results[counter,5]<-temp$coefficients[3]
                            results[counter,6]<-temp$additional[1]
                            results[counter,8]<-0
                            results[counter,9]<-temp$additional[4]
                        }
                        if(results[counter,8]>=80){
                            which_pos<-which(tumor[,2]>=start&tumor[,2]<=end)
                            which_inter<-intersect(which_chr,which_pos)
                            
                            keep_in[which_inter]<-T
                        }
                        counter<-counter+1
                    }
                }
            }
            dup[keep_in==F,2]<-NA
            progress$close()
            
            results<-results[order(results[,4],results[,1],results[,2]),]
            results20<-results[results[,8]>=80&((results[,4]=="del"&results[,6]>=(detectionThresholds[n,3]-0.05))|(results[,4]=="dup"&results[,6]>=(detectionThresholds[n,6]-0.05))),]
            if(sum(input$output_files=="Raw CNV calls")>0){
                write.table(results20,paste(input$output_folder,"/",samples_t[n,1],".CNVs_raw.txt",sep=""),row.names=F,sep="\t",quote=F)
            }

            if(sum(input$output_plots=="Raw CNV calls (all)")>0){
                output$text_plot1<-renderText({"Raw CNV calls (all)"})
                x.value_del<-data.frame(start=NA,end=NA,quality=NA,cells=NA)
                results_merged_del<-results[results[,4]=="del",]
                if(length(results_merged_del[,1])>0&&!is.na(results_merged_del[1,1])){
                    for(i in 1:length(results_merged_del[,1])){
                        if(results_merged_del[i,1]==1){
                            x.value_del[i,1]<-results_merged_del[i,2]
                            x.value_del[i,2]<-results_merged_del[i,3]
                        }
                        if(results_merged_del[i,1]>1){
                            x.value_del[i,1]<-genome[results_merged_del[i,1]-1]+results_merged_del[i,2]
                            x.value_del[i,2]<-genome[results_merged_del[i,1]-1]+results_merged_del[i,3]
                        }
                        x.value_del[i,3]<-log(abs(results_merged_del[i,8]))
                        x.value_del[i,4]<-results_merged_del[i,6]*(-1)
                    }
                }
                x.value_del[x.value_del[,3]==Inf|x.value_del[,3]==-Inf,3]<-NA
                x.value_del[is.na(x.value_del[,3]),3]<-min(max(x.value_del[,3],1,na.rm=T),100,na.rm=T)
                
                x.value_dup<-data.frame(start=NA,end=NA,quality=NA,cells=NA)
                results_merged_dup<-results[results[,4]=="dup",]
                if(length(results_merged_dup[,1])>0&&!is.na(results_merged_dup[1,1])){
                    for(i in 1:length(results_merged_dup[,1])){
                        if(results_merged_dup[i,1]==1){
                            x.value_dup[i,1]<-results_merged_dup[i,2]
                            x.value_dup[i,2]<-results_merged_dup[i,3]
                        }
                        if(results_merged_dup[i,1]>1){
                            x.value_dup[i,1]<-genome[results_merged_dup[i,1]-1]+results_merged_dup[i,2]
                            x.value_dup[i,2]<-genome[results_merged_dup[i,1]-1]+results_merged_dup[i,3]
                        }
                        x.value_dup[i,3]<-log(abs(results_merged_dup[i,8]))
                        x.value_dup[i,4]<-results_merged_dup[i,6]
                    }
                }
                x.value_dup[x.value_dup[,3]==Inf|x.value_dup[,3]==-Inf,3]<-NA
                x.value_dup[is.na(x.value_dup[,3]),3]<-min(max(x.value_dup[,3],1,na.rm=T),100,na.rm=T)
                
                png(paste0(input$output_folder,"/",samples_t[n,1],"_raw_all.png"),width=1800,height=800)
                plot(NULL,xlim=c(0,2881033286),ylim=c(-1.5,1.5),xaxt="n",yaxt="n",
                     xlab="Choromosome",ylab="",main=paste("Sample ",samples_t[n,1],sep=""))
                colfunc <- colorRampPalette(c("blue","red"))
                colpalette<-colfunc(max((max(x.value_del[,3])-min(x.value_del[,3]))*100+1,(max(x.value_dup[,3])-min(x.value_dup[,3]))*100+1))
                for(i in 1:length(x.value_del[,1])){
                    points(x.value_del[i,c(1,2)],c(x.value_del[i,4],x.value_del[i,4]),type="l",lwd=5,
                           col=colpalette[min(round(100*x.value_del[i,3])+1,length(colpalette))])
                }
                #colpalette<-colfunc((max(x.value_dup[,3])-min(x.value_dup[,3]))*100+1)
                for(i in 1:length(x.value_dup[,1])){
                    points(x.value_dup[i,c(1,2)],c(x.value_dup[i,4],x.value_dup[i,4]),type="l",lwd=5,
                           col=colpalette[min(round(100*x.value_dup[i,3])+1,length(colpalette))])
                }
                abline(h=0)
                abline(v=c(1,genome),col="grey")
                axis(1,at=c(1,genome),labels = NA)
                axis(1,at=c(c(c(1,genome[1:21])+genome)/2),labels = seq(1,22),col.ticks = NA)
                axis(2,at=c(-1,-0.5,0,0.5,1),labels = c("100%","50%","0%","50%","100%"))
                dev.off()
                output$plot1 <- renderImage({
                    list(src=paste0(input$output_folder,"/",samples_t[n,1],"_raw_all.png"),
                         height=400,
                         width=900)},
                    deleteFile = FALSE
                )  
            }
                
            if(length(results20[,1])>0){
                if(sum(input$output_plots=="Raw CNV calls (sig)")>0){
                    output$text_plot2<-renderText({"Raw CNV calls (sig)"})
                    #plot
                    no_dels<-T
                    x.value_del<-data.frame(start=NA,end=NA,quality=NA,cells=NA)
                    results_merged_del<-results20[results20[,4]=="del",]
                    if(length(results_merged_del[,1])>0&&!is.na(results_merged_del[1,1])){
                        for(i in 1:length(results_merged_del[,1])){
                            if(results_merged_del[i,1]==1){
                                x.value_del[i,1]<-results_merged_del[i,2]
                                x.value_del[i,2]<-results_merged_del[i,3]
                            }
                            if(results_merged_del[i,1]>1){
                                x.value_del[i,1]<-genome[results_merged_del[i,1]-1]+results_merged_del[i,2]
                                x.value_del[i,2]<-genome[results_merged_del[i,1]-1]+results_merged_del[i,3]
                            }
                            x.value_del[i,3]<-log(abs(results_merged_del[i,8]))
                            x.value_del[i,4]<-results_merged_del[i,6]*(-1)
                        }
                    }
                    if(length(x.value_del[,1])>1||!is.na(x.value_del[,1])){
                        x.value_del[x.value_del[,3]==Inf|x.value_del[,3]==-Inf,3]<-NA
                        x.value_del[is.na(x.value_del[,3]),3]<-min(max(x.value_del[,3],1,na.rm=T),100,na.rm=T)
                        no_dels<-F
                    }
                    no_dups<-T
                    x.value_dup<-data.frame(start=NA,end=NA,quality=NA,cells=NA)
                    results_merged_dup<-results20[results20[,4]=="dup",]
                    if(length(results_merged_dup[,1])>0&&!is.na(results_merged_dup[1,1])){
                        for(i in 1:length(results_merged_dup[,1])){
                            if(results_merged_dup[i,1]==1){
                                x.value_dup[i,1]<-results_merged_dup[i,2]
                                x.value_dup[i,2]<-results_merged_dup[i,3]
                            }
                            if(results_merged_dup[i,1]>1){
                                x.value_dup[i,1]<-genome[results_merged_dup[i,1]-1]+results_merged_dup[i,2]
                                x.value_dup[i,2]<-genome[results_merged_dup[i,1]-1]+results_merged_dup[i,3]
                            }
                            x.value_dup[i,3]<-log(abs(results_merged_dup[i,8]))
                            x.value_dup[i,4]<-results_merged_dup[i,6]
                        }
                    }
                    if(length(x.value_dup[,1])>1||!is.na(x.value_dup[,1])){
                        x.value_dup[x.value_dup[,3]==Inf|x.value_dup[,3]==-Inf,3]<-NA
                        x.value_dup[is.na(x.value_dup[,3]),3]<-min(max(x.value_dup[,3],1,na.rm=T),100,na.rm=T) 
                        no_dups<-F
                    }
                    
                    png(paste0(input$output_folder,"/",samples_t[n,1],"_raw_sig.png"),width=1800,height=800)
                    plot(NULL,xlim=c(0,2881033286),ylim=c(-1.5,1.5),xaxt="n",yaxt="n",
                         xlab="Choromosome",ylab="",main=paste("Sample ",samples_t[n,1],sep=""))
                    colfunc <- colorRampPalette(c("blue","red"))
                    if(no_dels==F&&no_dups==F){
                        colpalette<-colfunc(max((max(x.value_del[,3])-min(x.value_del[,3]))*100+1,(max(x.value_dup[,3])-min(x.value_dup[,3]))*100+1))
                    }
                    if(no_dels==F){
                        if(no_dups==T){
                            colpalette<-colfunc((max(x.value_del[,3])-min(x.value_del[,3]))*100+1)
                        }
                        for(i in 1:length(x.value_del[,1])){
                            points(x.value_del[i,c(1,2)],c(x.value_del[i,4],x.value_del[i,4]),type="l",lwd=5,
                                   col=colpalette[min(round(100*x.value_del[i,3])+1,length(colpalette))])
                        }
                    }
                    if(no_dups==F){
                        if(no_dels==T){
                            colpalette<-colfunc((max(x.value_dup[,3])-min(x.value_dup[,3]))*100+1)   
                        }
                        for(i in 1:length(x.value_dup[,1])){
                            points(x.value_dup[i,c(1,2)],c(x.value_dup[i,4],x.value_dup[i,4]),type="l",lwd=5,
                                   col=colpalette[min(round(100*x.value_dup[i,3])+1,length(colpalette))])
                        }   
                    }
                    abline(h=0)
                    abline(v=c(1,genome),col="grey")
                    axis(1,at=c(1,genome),labels = NA)
                    axis(1,at=c(c(c(1,genome[1:21])+genome)/2),labels = seq(1,22),col.ticks = NA)
                    axis(2,at=c(-1,-0.5,0,0.5,1),labels = c("100%","50%","0%","50%","100%"))
                    dev.off()
                    output$plot2 <- renderImage({
                        list(src=paste0(input$output_folder,"/",samples_t[n,1],"_raw_sig.png"),
                             height=400,
                             width=900)},
                        deleteFile = FALSE
                    )  
                }
                
                results20<-cbind(results20,windows=1)
                window2<-input$maxDist
                results_merged<-data.frame(chr=NA,start=NA,end=NA,variant=NA,p.value=NA,
                                           cells=NA,snps=NA,qual=NA,sd=NA,windows=NA,
                                           covTbelow=NA,covGbelow=NA)
                counter<-0
                
                if(length(results20[,1])>0&&!is.na(results20[1,1])){
                    shinyjs::html("text", paste0("&nbsp&nbsp&nbspMerging raw CNVs","<br>"), add = TRUE)
                    progress_sample$inc(1/4)
                    #merge: maximum 10 times
                    run_counter<-1
                    while(run_counter<=10){
                        #message(run_counter)
                        if(run_counter!=1){
                            results20<-results_merged
                            window2<-input$maxDist
                            results_merged<-data.frame(chr=NA,start=NA,end=NA,variant=NA,p.value=NA,
                                                       cells=NA,snps=NA,qual=NA,sd=NA,windows=NA,
                                                       covTbelow=NA,covGbelow=NA)
                            counter<-0
                        }
                        for(k in 1:(length(results20[,1]))){
                            if(!is.na(results20[k,5])){
                                #message(k)
                                counter<-counter+1
                                results_merged[counter,1]<-results20[k,1]
                                results_merged[counter,2]<-results20[k,2]
                                results_merged[counter,3]<-results20[k,3]
                                results_merged[counter,4]<-results20[k,4]
                                results_merged[counter,10]<-results20[k,10]
                                estimate<-results20[k,6]
                                estimate_low<-results20[k,6]-3*results20[k,9]
                                estimate_high<-results20[k,6]+3*results20[k,9]
                                
                                row.names(results20)<-seq(1,length(results20[,1]))
                                available_entries<-as.numeric(row.names(results20[results20[,1]==results_merged[counter,1]&!is.na(results20[,5]),]))
                                available_entries<-available_entries[available_entries>k]
                                for(i in available_entries){
                                    #overlapping; same variant -> test if similar
                                    if(results_merged[counter,3]>=results20[i,2]&&
                                       results_merged[counter,4]==results20[i,4]&&
                                       results_merged[counter,1]==results20[i,1]){
                                        
                                        chr<-results_merged[counter,1]
                                        tumor_1<-tumor[tumor[,1]==chr,]
                                        start<-results_merged[counter,2]
                                        end<-max(results20[i,3],results_merged[counter,3])
                                        
                                        if(results_merged[counter,4]=="del"){
                                            del_1<-del[tumor[,1]==chr,]
                                            del_g1<-del_g[tumor[,1]==chr,]
                                            del_oi<-del_1[tumor_1[,2]>=start&tumor_1[,2]<=end,]
                                            del_goi<-del_g1[tumor_1[,2]>=start&tumor_1[,2]<=end,]
                                            
                                            not_na<-intersect(which(!is.na(del_oi[,2])),which(!is.na(del_goi[,2])))
                                            del_oi<-del_oi[not_na,]
                                            del_goi<-del_goi[not_na,]
                                            temp<-wtd.t.test(x=del_oi[,2],y=del_goi[,2],
                                                             weight=del_oi[,6],weighty=del_goi[,6],
                                                             alternative = "greater")
                                        }
                                        if(results_merged[counter,4]=="dup"){
                                            dup_1<-dup[tumor[,1]==chr,]
                                            dup_g1<-dup_g[tumor[,1]==chr,]
                                            dup_oi<-dup_1[tumor_1[,2]>=start&tumor_1[,2]<=end,]
                                            dup_goi<-dup_g1[tumor_1[,2]>=start&tumor_1[,2]<=end,]
                                            
                                            not_na<-intersect(which(!is.na(dup_oi[,2])),which(!is.na(dup_goi[,2])))
                                            dup_oi<-dup_oi[not_na,]
                                            dup_goi<-dup_goi[not_na,]
                                            temp<-wtd.t.test(x=dup_oi[,2],y=dup_goi[,2],
                                                             weight=dup_oi[,6],weighty=dup_goi[,6],
                                                             alternative = "greater")
                                        }
                                        #inside -> merge
                                        
                                        
                                        if(sum(rowSums(cbind((temp$additional[1])<=estimate_high,
                                                             (temp$additional[1])>=estimate_low))>=1)==length(estimate_high)&&
                                           sum((temp$additional[1])<=(results20[i,6]+3*results20[i,9]))==1&&
                                           sum((temp$additional[1])>=(results20[i,6])-3*results20[i,9])==1){
                                            results_merged[counter,3]<-max(results20[i,3],results_merged[counter,3])
                                            results_merged[counter,10]<-sum(results_merged[counter,10],results20[i,10],na.rm=T)
                                            estimate<-c(estimate,results20[i,6])
                                            estimate_low<-temp$additional[1]-3*temp$additional[4]
                                            estimate_high<-temp$additional[1]+3*temp$additional[4]
                                            results20[i,5]<-NA
                                        }
                                    }
                                    #not overlapping, but same variant -> maybe merge
                                    if(results_merged[counter,3]<results20[i,2]&&
                                       results_merged[counter,4]==results20[i,4]&&
                                       results_merged[counter,1]==results20[i,1]){
                                        start_line<-(which(results_merged[counter,3]==results[,3]&results_merged[counter,1]==results[,1]&results_merged[counter,4]==results[,4])+1)[1]
                                        end_line<-(which(results20[i,2]==results[,2]&results20[i,1]==results[,1]&results20[i,4]==results[,4])-1)[1]
                                        if(end_line<start_line){
                                            #merging unproblematic; no regions in between
                                            chr<-results_merged[counter,1]
                                            tumor_1<-tumor[tumor[,1]==chr,]
                                            start<-results_merged[counter,2]
                                            end<-max(results20[i,3],results_merged[counter,3])
                                            
                                            if(results_merged[counter,4]=="del"){
                                                del_1<-del[tumor[,1]==chr,]
                                                del_g1<-del_g[tumor[,1]==chr,]
                                                del_oi<-del_1[tumor_1[,2]>=start&tumor_1[,2]<=end,]
                                                del_goi<-del_g1[tumor_1[,2]>=start&tumor_1[,2]<=end,]
                                                
                                                not_na<-intersect(which(!is.na(del_oi[,2])),which(!is.na(del_goi[,2])))
                                                del_oi<-del_oi[not_na,]
                                                del_goi<-del_goi[not_na,]
                                                temp<-wtd.t.test(x=del_oi[,2],y=del_goi[,2],
                                                                 weight=del_oi[,6],weighty=del_goi[,6],
                                                                 alternative = "greater")
                                            }
                                            if(results_merged[counter,4]=="dup"){
                                                dup_1<-dup[tumor[,1]==chr,]
                                                dup_g1<-dup_g[tumor[,1]==chr,]
                                                dup_oi<-dup_1[tumor_1[,2]>=start&tumor_1[,2]<=end,]
                                                dup_goi<-dup_g1[tumor_1[,2]>=start&tumor_1[,2]<=end,]
                                                
                                                not_na<-intersect(which(!is.na(dup_oi[,2])),which(!is.na(dup_goi[,2])))
                                                dup_oi<-dup_oi[not_na,]
                                                dup_goi<-dup_goi[not_na,]
                                                temp<-wtd.t.test(x=dup_oi[,2],y=dup_goi[,2],
                                                                 weight=dup_oi[,6],weighty=dup_goi[,6],
                                                                 alternative = "greater")
                                            }
                                            #inside -> merge
                                            if(sum(rowSums(cbind((temp$additional[1])<=estimate_high,
                                                                 (temp$additional[1])>=estimate_low))>=1)==length(estimate_high)&&
                                               sum((temp$additional[1])<=(results20[i,6]+3*results20[i,9]))==1&&
                                               sum((temp$additional[1])>=(results20[i,6])-3*results20[i,9])==1){
                                                results_merged[counter,3]<-max(results20[i,3],results_merged[counter,3])
                                                results_merged[counter,10]<-sum(results_merged[counter,10],results20[i,10],na.rm=T)
                                                estimate<-c(estimate,results20[i,6])
                                                estimate_low<-temp$additional[1]-3*temp$additional[4]
                                                estimate_high<-temp$additional[1]+3*temp$additional[4]
                                                results20[i,5]<-NA
                                            }
                                        }
                                        if(end_line>=start_line){
                                            #evaluate regions in between, then: decide if merge
                                            problematic<-results[start_line:end_line,]
                                            if(sum(problematic[,8]<=-80)==0&&(results20[i,2]-results_merged[counter,3])<window2){
                                                #no region significant for "less" and not too far away -> merge
                                                chr<-results_merged[counter,1]
                                                tumor_1<-tumor[tumor[,1]==chr,]
                                                start<-results_merged[counter,2]
                                                end<-max(results20[i,3],results_merged[counter,3])
                                                
                                                if(results_merged[counter,4]=="del"){
                                                    del_1<-del[tumor[,1]==chr,]
                                                    del_g1<-del_g[tumor[,1]==chr,]
                                                    del_oi<-del_1[tumor_1[,2]>=start&tumor_1[,2]<=end,]
                                                    del_goi<-del_g1[tumor_1[,2]>=start&tumor_1[,2]<=end,]
                                                    
                                                    not_na<-intersect(which(!is.na(del_oi[,2])),which(!is.na(del_goi[,2])))
                                                    del_oi<-del_oi[not_na,]
                                                    del_goi<-del_goi[not_na,]
                                                    temp<-wtd.t.test(x=del_oi[,2],y=del_goi[,2],
                                                                     weight=del_oi[,6],weighty=del_goi[,6],
                                                                     alternative = "greater")
                                                }
                                                if(results_merged[counter,4]=="dup"){
                                                    dup_1<-dup[tumor[,1]==chr,]
                                                    dup_g1<-dup_g[tumor[,1]==chr,]
                                                    dup_oi<-dup_1[tumor_1[,2]>=start&tumor_1[,2]<=end,]
                                                    dup_goi<-dup_g1[tumor_1[,2]>=start&tumor_1[,2]<=end,]
                                                    
                                                    not_na<-intersect(which(!is.na(dup_oi[,2])),which(!is.na(dup_goi[,2])))
                                                    dup_oi<-dup_oi[not_na,]
                                                    dup_goi<-dup_goi[not_na,]
                                                    temp<-wtd.t.test(x=dup_oi[,2],y=dup_goi[,2],
                                                                     weight=dup_oi[,6],weighty=dup_goi[,6],
                                                                     alternative = "greater")
                                                }
                                                #inside -> merge
                                                if(sum(rowSums(cbind((temp$additional[1])<=estimate_high,
                                                                     (temp$additional[1])>=estimate_low))>=1)==length(estimate_high)&&
                                                   sum((temp$additional[1])<=(results20[i,6]+3*results20[i,9]))==1&&
                                                   sum((temp$additional[1])>=(results20[i,6])-3*results20[i,9])==1){
                                                    results_merged[counter,3]<-max(results20[i,3],results_merged[counter,3])
                                                    results_merged[counter,10]<-sum(results_merged[counter,10],results20[i,10],na.rm=T)
                                                    estimate<-c(estimate,results20[i,6])
                                                    estimate_low<-temp$additional[1]-3*temp$additional[4]
                                                    estimate_high<-temp$additional[1]+3*temp$additional[4]
                                                    results20[i,5]<-NA
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        if(length(grep("del",results_merged[,4]))>0){
                            for(j in grep("del",results_merged[,4])){
                                chr<-results_merged[j,1]
                                tumor_1<-tumor[tumor[,1]==chr,]
                                del_1<-del[tumor[,1]==chr,]
                                del_g1<-del_g[tumor[,1]==chr,]
                                
                                start<-results_merged[j,2]
                                end<-results_merged[j,3]
                                del_oi<-del_1[tumor_1[,2]>=start&tumor_1[,2]<=end,]
                                del_goi<-del_g1[tumor_1[,2]>=start&tumor_1[,2]<=end,]
                                
                                not_na<-intersect(which(!is.na(del_oi[,2])),which(!is.na(del_goi[,2])))
                                del_oi<-del_oi[not_na,]
                                del_goi<-del_goi[not_na,]
                                temp<-wtd.t.test(x=del_oi[,2],y=del_goi[,2],
                                                 weight=del_oi[,6],weighty=del_goi[,6],
                                                 alternative = "greater")
                                results_merged[j,7]<-length(del_goi[,1])
                                
                                results_merged[j,5]<-temp$coefficients[3]
                                results_merged[j,6]<-temp$additional[1]
                                results_merged[j,8]<-temp$coefficients[3]**(-1)
                                results_merged[j,9]<-temp$additional[4]
                                results_merged[j,11]<-sum(del_oi[,4],na.rm=T)
                                results_merged[j,12]<-sum(del_goi[,4],na.rm=T)
                            }
                        }
                        
                        if(length(grep("dup",results_merged[,4]))>0){
                            for(j in grep("dup",results_merged[,4])){
                                chr<-results_merged[j,1]
                                tumor_1<-tumor[tumor[,1]==chr,]
                                dup_1<-dup[tumor[,1]==chr,]
                                dup_g1<-dup_g[tumor[,1]==chr,]
                                
                                start<-results_merged[j,2]
                                end<-results_merged[j,3]
                                dup_oi<-dup_1[tumor_1[,2]>=start&tumor_1[,2]<=end,]
                                dup_goi<-dup_g1[tumor_1[,2]>=start&tumor_1[,2]<=end,]
                                
                                not_na<-intersect(which(!is.na(dup_oi[,2])),which(!is.na(dup_goi[,2])))
                                dup_oi<-dup_oi[not_na,]
                                dup_goi<-dup_goi[not_na,]
                                temp<-wtd.t.test(x=dup_oi[,2],y=dup_goi[,2],
                                                 weight=dup_oi[,6],weighty=dup_goi[,6],
                                                 alternative = "greater")
                                results_merged[j,7]<-sum(!is.na(dup_oi[,1]))
                                
                                results_merged[j,5]<-temp$coefficients[3]
                                results_merged[j,6]<-temp$additional[1]
                                results_merged[j,8]<-temp$coefficients[3]**(-1)
                                results_merged[j,9]<-temp$additional[4]
                                results_merged[j,11]<-sum(dup_oi[,4],na.rm=T)
                                results_merged[j,12]<-sum(dup_goi[!is.na(dup_oi[,1]),4],na.rm=T)
                            }
                        }
                        
                        if(length(results_merged[,1])==length(results20[,1])){
                            run_counter<-11
                        }
                        if(length(results_merged[,1])!=length(results20[,1])){
                            run_counter<-run_counter+1
                        }
                    }
                    
                    results_merged<-results_merged[results_merged[,8]>=80,]
                    results_merged<-results_merged[results_merged[,6]<1.5,]

                    if(length(results_merged[,1])>0&&!is.na(results_merged[1,1])){
                        results_merged[,5]<-round(results_merged[,5],digits = 4)
                        results_merged[,6]<-round(results_merged[,6],digits = 8)
                        results_merged[,8]<-round(results_merged[,8],digits = 2)
                        results_merged[,9]<-round(results_merged[,9],digits = 2)
                        
                        results_merged<-cbind(results_merged,logQual=log(results_merged[,8]))
                        results_merged[,13]<-round(results_merged[,13],digits = 2)
                        results_merged<-cbind(results_merged,covIndicator=NA)
                        for(i in 1:length(results_merged[,1])){
                            cov_indicator_1<-binom.test(x=results_merged[i,11],n=results_merged[i,7],
                                                        p=results_merged[i,12]/results_merged[i,7])
                            if(round(results_merged[i,12]/results_merged[i,7]-cov_indicator_1$estimate,2)>0){
                                results_merged[i,14]<-paste0("+",round(results_merged[i,12]/results_merged[i,7]-cov_indicator_1$estimate,2)," [",
                                                             round(results_merged[i,12]/results_merged[i,7]-cov_indicator_1$conf.int[1],2),";",
                                                             round(results_merged[i,12]/results_merged[i,7]-cov_indicator_1$conf.int[2],2),"]")          
                            }
                            if(round(results_merged[i,12]/results_merged[i,7]-cov_indicator_1$estimate,2)<=0){
                                results_merged[i,14]<-paste0(round(results_merged[i,12]/results_merged[i,7]-cov_indicator_1$estimate,2)," [",
                                                             round(results_merged[i,12]/results_merged[i,7]-cov_indicator_1$conf.int[1],2),";",
                                                             round(results_merged[i,12]/results_merged[i,7]-cov_indicator_1$conf.int[2],2),"]")          
                            }
                        }
                        
                        results_merged<-cbind(Sample=samples_t[n,1],results_merged,stringsAsFactors=F)
                        results_merged_report<-results_merged[,c(1:7,10,8,11,14,15)]
                        names(results_merged_report)<-c("Sample","chr","start","end","variant","p value","estimated CF",
                                                        "SD","SNPs","Windows","Quality","Cov-indicator")
                        
                        if(sum(input$output_files=="Merged CNV calls")>0){
                            write.table(results_merged_report,paste0(input$output_folder,"/",samples_t[n,1],".CNVs_merged.txt"),row.names=F,sep="\t",quote=F)   
                        }
                        
                        if(sum(input$output_plots=="Merged CNV calls")>0){
                            output$text_plot3<-renderText({"Merged CNV calls"})
                            #plot
                            no_dels<-T
                            x.value_del<-data.frame(start=NA,end=NA,quality=NA,cells=NA)
                            results_merged_del<-results_merged[results_merged[,5]=="del",c(2:15)]
                            if(length(results_merged_del[,1])>0&&!is.na(results_merged_del[1,1])){
                                for(i in 1:length(results_merged_del[,1])){
                                    if(results_merged_del[i,1]==1){
                                        x.value_del[i,1]<-results_merged_del[i,2]
                                        x.value_del[i,2]<-results_merged_del[i,3]
                                    }
                                    if(results_merged_del[i,1]>1){
                                        x.value_del[i,1]<-genome[results_merged_del[i,1]-1]+results_merged_del[i,2]
                                        x.value_del[i,2]<-genome[results_merged_del[i,1]-1]+results_merged_del[i,3]
                                    }
                                    x.value_del[i,3]<-results_merged_del[i,13]
                                    x.value_del[i,4]<-results_merged_del[i,6]*(-1)
                                }
                            }
                            if(length(x.value_del[,1])>1||!is.na(x.value_del[,1])){
                                x.value_del[x.value_del[,3]==Inf|x.value_del[,3]==-Inf,3]<-NA
                                x.value_del[is.na(x.value_del[,3]),3]<-min(max(x.value_del[,3],1,na.rm=T),100,na.rm=T)
                                no_dels<-F
                            }
                            no_dups<-T
                            x.value_dup<-data.frame(start=NA,end=NA,quality=NA,cells=NA)
                            results_merged_dup<-results_merged[results_merged[,5]=="dup",c(2:15)]
                            if(length(results_merged_dup[,1])>0&&!is.na(results_merged_dup[1,1])){
                                for(i in 1:length(results_merged_dup[,1])){
                                    if(results_merged_dup[i,1]==1){
                                        x.value_dup[i,1]<-results_merged_dup[i,2]
                                        x.value_dup[i,2]<-results_merged_dup[i,3]
                                    }
                                    if(results_merged_dup[i,1]>1){
                                        x.value_dup[i,1]<-genome[results_merged_dup[i,1]-1]+results_merged_dup[i,2]
                                        x.value_dup[i,2]<-genome[results_merged_dup[i,1]-1]+results_merged_dup[i,3]
                                    }
                                    x.value_dup[i,3]<-results_merged_dup[i,13]
                                    x.value_dup[i,4]<-results_merged_dup[i,6]
                                }
                            }
                            if(length(x.value_dup[,1])>1||!is.na(x.value_dup[,1])){
                                x.value_dup[x.value_dup[,3]==Inf|x.value_dup[,3]==-Inf,3]<-NA
                                x.value_dup[is.na(x.value_dup[,3]),3]<-min(max(x.value_dup[,3],1,na.rm=T),100,na.rm=T) 
                                no_dups<-F
                            }
                            png(paste0(input$output_folder,"/",samples_t[n,1],"_merged.png"),width=1800,height=800)
                            plot(NULL,xlim=c(0,2881033286),ylim=c(-1.5,1.5),xaxt="n",yaxt="n",
                                 xlab="Choromosome",ylab="",main=paste("Sample ",samples_t[n,1],sep=""))
                            colfunc <- colorRampPalette(c("blue","red"))
                            if(no_dels==F&&no_dups==F){
                                colpalette<-colfunc(max((max(x.value_del[,3])-min(x.value_del[,3]))*100+1,(max(x.value_dup[,3])-min(x.value_dup[,3]))*100+1))
                            }
                            if(no_dels==F){
                                if(no_dups==T){
                                    colpalette<-colfunc((max(x.value_del[,3])-min(x.value_del[,3]))*100+1)
                                }
                                for(i in 1:length(x.value_del[,1])){
                                    points(x.value_del[i,c(1,2)],c(x.value_del[i,4],x.value_del[i,4]),type="l",lwd=5,
                                           col=colpalette[min(round(100*x.value_del[i,3])+1,length(colpalette))])
                                }            
                            }
                            if(no_dups==F){
                                if(no_dels==T){
                                    colpalette<-colfunc((max(x.value_dup[,3])-min(x.value_dup[,3]))*100+1)
                                }
                                for(i in 1:length(x.value_dup[,1])){
                                    points(x.value_dup[i,c(1,2)],c(x.value_dup[i,4],x.value_dup[i,4]),type="l",lwd=5,
                                           col=colpalette[min(round(100*x.value_dup[i,3])+1,length(colpalette))])
                                }   
                            }
                            abline(h=0)
                            abline(v=c(1,genome),col="grey")
                            axis(1,at=c(1,genome),labels = NA)
                            axis(1,at=c(c(c(1,genome[1:21])+genome)/2),labels = seq(1,22),col.ticks = NA)
                            axis(2,at=c(-1,-0.5,0,0.5,1),labels = c("100%","50%","0%","50%","100%"))
                            dev.off()
                            output$plot3 <- renderImage({
                                list(src=paste0(input$output_folder,"/",samples_t[n,1],"_merged.png"),
                                     height=400,
                                     width=900)},
                                deleteFile = FALSE
                            )  
                        }
                        
                        if(input$final_filter=="No"){
                            output$table_cnvs <- renderDataTable(datatable(results_merged_report))
                        }
                        if(input$final_filter=="Yes"){
                            shinyjs::html("text", paste0("&nbsp&nbsp&nbspFilter merged CNVs","<br>"), add = TRUE)
                            progress_sample$inc(1/4)
                            results_final<-results_merged[results_merged[,14]>=input$quality_filter,]
                            results_final_report<-results_final[,c(1:7,10,8,11,14,15)]
                            names(results_final_report)<-c("Sample","chr","start","end","variant","p value","estimated CF",
                                                           "SD","SNPs","Windows","Quality","Cov-indicator")
                            output$table_cnvs <- renderDataTable(datatable(results_final_report))
                            if(sum(input$output_files=="Filtered CNV calls")>0){
                                write.table(results_final_report,paste0(input$output_folder,"/",samples_t[n,1],".CNVs_filtered.txt"),row.names=F,sep="\t",quote=F)   
                            }
                            
                            if(sum(input$output_plots=="Filtered CNV calls")>0){
                                output$text_plot4<-renderText({"Filtered CNV calls"})
                                #plot
                                no_dels<-T
                                x.value_del<-data.frame(start=NA,end=NA,quality=NA,cells=NA)
                                results_final_del<-results_final[results_final[,5]=="del",c(2:15)]
                                if(length(results_final_del[,1])>0&&!is.na(results_final_del[1,1])){
                                    for(i in 1:length(results_final_del[,1])){
                                        if(results_final_del[i,1]==1){
                                            x.value_del[i,1]<-results_final_del[i,2]
                                            x.value_del[i,2]<-results_final_del[i,3]
                                        }
                                        if(results_final_del[i,1]>1){
                                            x.value_del[i,1]<-genome[results_final_del[i,1]-1]+results_final_del[i,2]
                                            x.value_del[i,2]<-genome[results_final_del[i,1]-1]+results_final_del[i,3]
                                        }
                                        x.value_del[i,3]<-results_final_del[i,13]
                                        x.value_del[i,4]<-results_final_del[i,6]*(-1)
                                    }
                                }
                                if(length(x.value_del[,1])>1||!is.na(x.value_del[,1])){
                                    x.value_del[x.value_del[,3]==Inf|x.value_del[,3]==-Inf,3]<-NA
                                    x.value_del[is.na(x.value_del[,3]),3]<-min(max(x.value_del[,3],1,na.rm=T),100,na.rm=T)
                                    no_dels<-F
                                }
                                no_dups<-T
                                x.value_dup<-data.frame(start=NA,end=NA,quality=NA,cells=NA)
                                results_final_dup<-results_final[results_final[,5]=="dup",c(2:15)]
                                if(length(results_final_dup[,1])>0&&!is.na(results_final_dup[1,1])){
                                    for(i in 1:length(results_final_dup[,1])){
                                        if(results_final_dup[i,1]==1){
                                            x.value_dup[i,1]<-results_final_dup[i,2]
                                            x.value_dup[i,2]<-results_final_dup[i,3]
                                        }
                                        if(results_final_dup[i,1]>1){
                                            x.value_dup[i,1]<-genome[results_final_dup[i,1]-1]+results_final_dup[i,2]
                                            x.value_dup[i,2]<-genome[results_final_dup[i,1]-1]+results_final_dup[i,3]
                                        }
                                        x.value_dup[i,3]<-results_final_dup[i,13]
                                        x.value_dup[i,4]<-results_final_dup[i,6]
                                    }
                                }
                                if(length(x.value_dup[,1])>1||!is.na(x.value_dup[,1])){
                                    x.value_dup[x.value_dup[,3]==Inf|x.value_dup[,3]==-Inf,3]<-NA
                                    x.value_dup[is.na(x.value_dup[,3]),3]<-min(max(x.value_dup[,3],1,na.rm=T),100,na.rm=T)
                                    no_dups<-F
                                }
                                
                                png(paste0(input$output_folder,"/",samples_t[n,1],"_filtered.png"),width=1800,height=800)
                                plot(NULL,xlim=c(0,2881033286),ylim=c(-1.5,1.5),xaxt="n",yaxt="n",
                                     xlab="Choromosome",ylab="",main=paste("Sample ",samples_t[n,1],sep=""))
                                colfunc <- colorRampPalette(c("blue","red"))
                                if(no_dels==F&&no_dups==F){
                                    colpalette<-colfunc(max((max(x.value_del[,3])-min(x.value_del[,3]))*100+1,(max(x.value_dup[,3])-min(x.value_dup[,3]))*100+1))
                                }
                                if(no_dels==F){
                                    if(no_dups==T){
                                        colpalette<-colfunc((max(x.value_del[,3])-min(x.value_del[,3]))*100+1)
                                    }
                                    for(i in 1:length(x.value_del[,1])){
                                        points(x.value_del[i,c(1,2)],c(x.value_del[i,4],x.value_del[i,4]),type="l",lwd=5,
                                               col=colpalette[min(round(100*x.value_del[i,3])+1,length(colpalette))])
                                    }            
                                }
                                if(no_dups==F){
                                    if(no_dels==T){
                                        colpalette<-colfunc((max(x.value_dup[,3])-min(x.value_dup[,3]))*100+1)
                                    }
                                    for(i in 1:length(x.value_dup[,1])){
                                        points(x.value_dup[i,c(1,2)],c(x.value_dup[i,4],x.value_dup[i,4]),type="l",lwd=5,
                                               col=colpalette[min(round(100*x.value_dup[i,3])+1,length(colpalette))])
                                    }   
                                }
                                abline(h=0)
                                abline(v=c(1,genome),col="grey")
                                axis(1,at=c(1,genome),labels = NA)
                                axis(1,at=c(c(c(1,genome[1:21])+genome)/2),labels = seq(1,22),col.ticks = NA)
                                axis(2,at=c(-1,-0.5,0,0.5,1),labels = c("100%","50%","0%","50%","100%"))
                                dev.off()
                                output$plot4 <- renderImage({
                                    list(src=paste0(input$output_folder,"/",samples_t[n,1],"_filtered.png"),
                                         height=400,
                                         width=900)},
                                    deleteFile = FALSE
                                )  
                            }
                            
                        }
                    }
                }
            }
           
            
            if(length(results20[,1])==0){
                shinyjs::html("text", paste0("<br>","&nbsp&nbsp&nbspNo CNVs detected","<br>"), add = TRUE)
            }
            progress_sample$close()
            write.table(erwartungswerte,paste0(input$output_folder,"/Erwartungswerte.txt"),
                        sep="\t",row.names = F,quote=F)
        }
        updateRadioButtons(session,"select_samples",choices=samples_t[,1],
                           selected = samples_t[length(samples_t[,1]),1],inline=T)
        
        write.table(detectionThresholds,paste0(input$output_folder,"/DetectionThresholds.txt"),
                    sep="\t",row.names = F,quote=F)
        
        write.table(erwartungswerte,paste0(input$output_folder,"/Erwartungswerte.txt"),
                    sep="\t",row.names = F,quote=F)
        
        shinyjs::html("text", paste0("<br>","Analysis complete","<br>"), add = TRUE)
    })
    
    #Display results
    observeEvent(input$do2,{
        if(input$select_samples[1]==""){
            shinyjs::html("text", paste0("<br>","Please perform analysis first","<br>"), add = FALSE)
        }
        if(!is.na(input$select_samples)){
            output$sample<-renderText({paste0("Sample ",input$select_samples)})
            output$sample2<-renderText({paste0("Sample ",input$select_samples)})
            
            output$plot1<-NULL
            output$plot2<-NULL
            output$plot3<-NULL
            output$plot4<-NULL
            output$text_plot1<-renderText({""})
            output$text_plot2<-renderText({""})
            output$text_plot3<-renderText({""})
            output$text_plot4<-renderText({""})
            if(sum(input$output_plots2=="Raw CNV calls (all)")>0&&
               file.exists(paste0(input$output_folder,"/",input$select_samples,"_raw_all.png"))){
                output$text_plot1<-renderText({"Raw CNV calls (all)"})
                output$plot1 <- renderImage({
                    list(src=paste0(input$output_folder,"/",input$select_samples,"_raw_all.png"),
                         height=400,
                         width=900)},
                    deleteFile = FALSE
                )  
            }
            if(sum(input$output_plots2=="Raw CNV calls (sig)")>0&&
               file.exists(paste0(input$output_folder,"/",input$select_samples,"_raw_sig.png"))){
                if(sum(input$output_plots2=="Raw CNV calls (all)")>0){
                    output$text_plot2<-renderText({"Raw CNV calls (sig)"})
                    output$plot2 <- renderImage({
                        list(src=paste0(input$output_folder,"/",input$select_samples,"_raw_sig.png"),
                             height=400,
                             width=900)},
                        deleteFile = FALSE
                    )  
                }
                if(sum(input$output_plots2=="Raw CNV calls (all)")==0){
                    output$text_plot1<-renderText({"Raw CNV calls (sig)"})
                    output$plot1 <- renderImage({
                        list(src=paste0(input$output_folder,"/",input$select_samples,"_raw_sig.png"),
                             height=400,
                             width=900)},
                        deleteFile = FALSE
                    )  
                }

            }
            if(sum(input$output_plots2=="Merged CNV calls")>0&&
               file.exists(paste0(input$output_folder,"/",input$select_samples,"_merged.png"))){
                if(sum(input$output_plots2=="Raw CNV calls (all)")>0){
                    if(sum(input$output_plots2=="Raw CNV calls (sig)")>0){
                        output$text_plot3<-renderText({"Merged CNV calls"})
                        output$plot3 <- renderImage({
                            list(src=paste0(input$output_folder,"/",input$select_samples,"_merged.png"),
                                 height=400,
                                 width=900)},
                            deleteFile = FALSE
                        ) 
                    }
                    if(sum(input$output_plots2=="Raw CNV calls (sig)")==0){
                        output$text_plot2<-renderText({"Merged CNV calls"})
                        output$plot2 <- renderImage({
                            list(src=paste0(input$output_folder,"/",input$select_samples,"_merged.png"),
                                 height=400,
                                 width=900)},
                            deleteFile = FALSE
                        ) 
                    }
                }
                if(sum(input$output_plots2=="Raw CNV calls (all)")==0){
                    if(sum(input$output_plots2=="Raw CNV calls (sig)")>0){
                        output$text_plot2<-renderText({"Merged CNV calls"})
                        output$plot2 <- renderImage({
                            list(src=paste0(input$output_folder,"/",input$select_samples,"_merged.png"),
                                 height=400,
                                 width=900)},
                            deleteFile = FALSE
                        ) 
                    }
                    if(sum(input$output_plots2=="Raw CNV calls (sig)")==0){
                        output$text_plot1<-renderText({"Merged CNV calls"})
                        output$plot1 <- renderImage({
                            list(src=paste0(input$output_folder,"/",input$select_samples,"_merged.png"),
                                 height=400,
                                 width=900)},
                            deleteFile = FALSE
                        ) 
                    }
                }
 
            }
            if(sum(input$output_plots2=="Filtered CNV calls")>0&&
               file.exists(paste0(input$output_folder,"/",input$select_samples,"_filtered.png"))){
                if(sum(input$output_plots2=="Raw CNV calls (all)")>0){
                    if(sum(input$output_plots2=="Raw CNV calls (sig)")>0){
                        if(sum(input$output_plots2=="Merged CNV calls")>0){
                            output$text_plot4<-renderText({"Filtered CNV calls"})
                            output$plot4 <- renderImage({
                                list(src=paste0(input$output_folder,"/",input$select_samples,"_filtered.png"),
                                     height=400,
                                     width=900)},
                                deleteFile = FALSE
                            ) 
                        }
                        if(sum(input$output_plots2=="Merged CNV calls")==0){
                            output$text_plot3<-renderText({"Filtered CNV calls"})
                            output$plot3 <- renderImage({
                                list(src=paste0(input$output_folder,"/",input$select_samples,"_filtered.png"),
                                     height=400,
                                     width=900)},
                                deleteFile = FALSE
                            ) 
                        }
                    }
                    if(sum(input$output_plots2=="Raw CNV calls (sig)")==0){
                        if(sum(input$output_plots2=="Merged CNV calls")>0){
                            output$text_plot3<-renderText({"Filtered CNV calls"})
                            output$plot3 <- renderImage({
                                list(src=paste0(input$output_folder,"/",input$select_samples,"_filtered.png"),
                                     height=400,
                                     width=900)},
                                deleteFile = FALSE
                            ) 
                        }
                        if(sum(input$output_plots2=="Merged CNV calls")==0){
                            output$text_plot2<-renderText({"Filtered CNV calls"})
                            output$plot2 <- renderImage({
                                list(src=paste0(input$output_folder,"/",input$select_samples,"_filtered.png"),
                                     height=400,
                                     width=900)},
                                deleteFile = FALSE
                            ) 
                        }
                    }
                }
                if(sum(input$output_plots2=="Raw CNV calls (all)")==0){
                    if(sum(input$output_plots2=="Raw CNV calls (sig)")>0){
                        if(sum(input$output_plots2=="Merged CNV calls")>0){
                            output$text_plot3<-renderText({"Filtered CNV calls"})
                            output$plot3 <- renderImage({
                                list(src=paste0(input$output_folder,"/",input$select_samples,"_filtered.png"),
                                     height=400,
                                     width=900)},
                                deleteFile = FALSE
                            ) 
                        }
                        if(sum(input$output_plots2=="Merged CNV calls")==0){
                            output$text_plot2<-renderText({"Filtered CNV calls"})
                            output$plot2 <- renderImage({
                                list(src=paste0(input$output_folder,"/",input$select_samples,"_filtered.png"),
                                     height=400,
                                     width=900)},
                                deleteFile = FALSE
                            ) 
                        }
                    }
                    if(sum(input$output_plots2=="Raw CNV calls (sig)")==0){
                        if(sum(input$output_plots2=="Merged CNV calls")>0){
                            output$text_plot2<-renderText({"Filtered CNV calls"})
                            output$plot2 <- renderImage({
                                list(src=paste0(input$output_folder,"/",input$select_samples,"_filtered.png"),
                                     height=400,
                                     width=900)},
                                deleteFile = FALSE
                            ) 
                        }
                        if(sum(input$output_plots2=="Merged CNV calls")==0){
                            output$text_plot1<-renderText({"Filtered CNV calls"})
                            output$plot1 <- renderImage({
                                list(src=paste0(input$output_folder,"/",input$select_samples,"_filtered.png"),
                                     height=400,
                                     width=900)},
                                deleteFile = FALSE
                            ) 
                        }
                    }
                }
 
            }
            
            if(input$output_files2=="Merged CNV calls"&&
               file.exists(paste0(input$output_folder,"/",input$select_samples,".CNVs_merged.txt"))){
                calls<-read.table(paste0(input$output_folder,"/",input$select_samples,".CNVs_merged.txt"),
                                  header=T,quote = "",sep="\t",stringsAsFactors = F)
                #calls[,13]<-as.character(calls[,12])
                output$table_cnvs <- renderDataTable(datatable(calls))
            }
            if(input$output_files2=="Filtered CNV calls"&&
               file.exists(paste0(input$output_folder,"/",input$select_samples,".CNVs_filtered.txt"))){
                calls<-read.table(paste0(input$output_folder,"/",input$select_samples,".CNVs_filtered.txt"),
                                  header=T,quote = "",sep="\t",stringsAsFactors = F)
                #calls[,13]<-as.character(calls[,13])
                output$table_cnvs <- renderDataTable(datatable(calls))
            }
        }

    })
}
)
