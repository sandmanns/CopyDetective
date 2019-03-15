shinyServer(function(input, output, session) {
    
    output$simulationUI<-renderUI({conditionalPanel(
        condition="input.simulation=='Yes'",
        fileInput('detectionThresholdsFile',label = "Upload detection thresholds file")
    )})
    output$simulationUI2<-renderUI({conditionalPanel(
        condition="input.simulation=='No'",
        numericInput('simulations_nr',"Number of simulations to perform [10, 99999]",min=10,max=99999,value=500),
        hr(),
        h4("A: Evaluate Tumor Ratios:"),
        sliderInput('simulate_ratios',"Ratios to consider [%]",min=1,max=100,value=c(5,100)),
        numericInput('simulate_steps',"In steps of...",min=1,max=100,value=5),
        hr(),
        h4("B: Evaluate Window Sizes:"),
        numericInput('simulate_window_start',"Minimum window size to consider [100000, 59128983]",min=100000,max=59128983,value=1000000),
        numericInput('simulate_window_steps',"In steps of... [100000, 59128983]",min=100000,max=59128983,value=1000000),
        hr(),
        radioButtons('simulation_optimization',label="Select strategy for detection threshold optimization",
                     choices=c("Compromize (ratio and window size)","Force (ratio)","Force (window size)"),
                     selected = "Compromize (ratio and window size)",inline = T)
    )})
    output$simulationUI3<-renderUI({conditionalPanel(
        condition="input.simulation_optimization=='Force (ratio)'",
        numericInput('simulations_force_ratio',"Minimum ratio to consider [%]",min=1,max=100,value=10)
    )})
    output$simulationUI4<-renderUI({conditionalPanel(
        condition="input.simulation_optimization=='Force (window size)'",
        numericInput('simulations_force_window',"Minimum window size to consider [bp]",min=100000,max=59128983,value=1000000)
    )})
    
    output$filtration_thresholdUI<-renderUI({conditionalPanel(
        condition="input.final_filter=='Yes'",
        numericInput('quality_filter','Quality Threshold',value=6.00,min=0,max=50,step = 0.01)
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
        dir2<-input$input_folder
        dir_out<-input$output_folder
        
        samples_temp<-input$sampleFile
        if(is.null(samples_temp)){
            shinyjs::html("text", paste0("No sample names file defined","<br>"), add = TRUE) 
            return()
        }
        samples<-read.table(samples_temp$datapath,header=F,quote = "",sep="\t",stringsAsFactors = F)
        if(sum(is.na(samples[,1]))!=0||sum(is.na(samples[,2]))!=0){
            shinyjs::html("text", paste0("At least one sample without matching 2nd sample","<br>"), add = TRUE) 
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
        
        
        detectionThresholds<-data.frame(Sample=samples_t[,1],Window_del=NA,Ratio_del=NA,
                                        Window_dup=NA,Ratio_dup=NA)

        for(n in 1:length(samples_t[,1])){
            shinyjs::html("text", paste0("<br>Starting analysis with CopyDetective...<br><br>"), add = FALSE)
            shinyjs::html("text", paste0("Analyzing sample ",samples_t[n,1],"<br>"), add = TRUE)
            progress_sample <- shiny::Progress$new()
            progress_sample$set(message = paste0("Analyzing sample ",samples_t[n,1]), value = 0)
            output$sample<-renderText({paste0("Sample ",samples_t[n,1])})
            output$sample2<-renderText({paste0("Sample ",samples_t[n,1])})
            
            germline<-read.table(paste0(input$input_folder,samples_g[n,1],".txt"),header = T,sep="\t",stringsAsFactors = F)
            tumor<-read.table(paste0(input$input_folder,samples_t[n,1],".txt"),header = T,sep="\t",stringsAsFactors = F)
            germline<-cbind(germline[,c(1:4)],Sample=samples_g[n,1],germline[,c(5:8)],stringsAsFactors = F)
            tumor<-cbind(tumor[,c(1:4)],Sample=samples_t[n,1],tumor[,c(5:8)],stringsAsFactors = F)
            
            tumor<-tumor[order(tumor[,1],tumor[,2]),]
            germline<-germline[order(germline[,1],germline[,2]),]
            
            #Preparations: determine the cellular ratios and the confidence intervals
            #for tumor
            del<-data.frame(upper=NA,mean=NA,lower=NA)
            dup<-data.frame(upper=NA,mean=NA,lower=NA)
            for(i in 1:length(tumor[,1])){
                if(!is.na(tumor[i,9])){
                    lower<-binom.test(x=tumor[i,7],tumor[i,8],tumor[i,9],conf.level = 0.95)$conf.int[1]
                    upper<-binom.test(x=tumor[i,7],tumor[i,8],tumor[i,9],conf.level = 0.95)$conf.int[2]
                    if(tumor[i,9]<=0.5){
                        del[i,1]<-(2*lower-1)/(lower-1)
                        del[i,2]<-(2*tumor[i,9]-1)/(tumor[i,9]-1)
                        del[i,3]<-(2*upper-1)/(upper-1)
                        
                        if(upper>=1/3&&tumor[i,9]!=0){
                            dup[i,1]<-(2*lower-1)/(-1*lower)
                            dup[i,2]<-(2*tumor[i,9]-1)/(-1*tumor[i,9])
                            dup[i,3]<-(2*upper-1)/(-1*upper)
                        }
                    }
                    if(tumor[i,9]>0.5){
                        del[i,3]<-(2*lower-1)/(lower)
                        del[i,2]<-(2*tumor[i,9]-1)/(tumor[i,9])
                        del[i,1]<-(2*upper-1)/(upper)
                        
                        if(lower<=2/3&&tumor[i,9]!=1){
                            dup[i,3]<-(2*lower-1)/(1-lower)
                            dup[i,2]<-(2*tumor[i,9]-1)/(1-tumor[i,9])
                            dup[i,1]<-(2*upper-1)/(1-upper)
                        }
                    }
                }
            }
            tumor[is.na(tumor[,8]),8]<-1
            
            #for germline
            del_g<-data.frame(upper=NA,mean=NA,lower=NA)
            dup_g<-data.frame(upper=NA,mean=NA,lower=NA)
            for(i in 1:length(tumor[,1])){
                if(!is.na(germline[i,9])){
                    lower<-binom.test(x=germline[i,7],germline[i,8],germline[i,9],conf.level = 0.95)$conf.int[1]
                    upper<-binom.test(x=germline[i,7],germline[i,8],germline[i,9],conf.level = 0.95)$conf.int[2]
                    if(germline[i,9]<=0.5){
                        del_g[i,1]<-(2*lower-1)/(lower-1)
                        del_g[i,2]<-(2*germline[i,9]-1)/(germline[i,9]-1)
                        del_g[i,3]<-(2*upper-1)/(upper-1)
                        
                        if(upper>=1/3&&germline[i,9]!=0){
                            dup_g[i,1]<-(2*lower-1)/(-1*lower)
                            dup_g[i,2]<-(2*germline[i,9]-1)/(-1*germline[i,9])
                            dup_g[i,3]<-(2*upper-1)/(-1*upper)
                        }
                    }
                    if(germline[i,9]>0.5){
                        del_g[i,3]<-(2*lower-1)/(lower)
                        del_g[i,2]<-(2*germline[i,9]-1)/(germline[i,9])
                        del_g[i,1]<-(2*upper-1)/(upper)
                        
                        if(lower<=2/3&&germline[i,9]!=1){
                            dup_g[i,3]<-(2*lower-1)/(1-lower)
                            dup_g[i,2]<-(2*germline[i,9]-1)/(1-germline[i,9])
                            dup_g[i,1]<-(2*upper-1)/(1-upper)
                        }
                    }
                }
            }
            
            progress_sample$inc(1/4)
            if(input$simulation=="No"){
                shinyjs::html("text", paste0("<br>","&nbsp&nbsp&nbspDetermining thresholds for CNV detection","<br>"), add = TRUE)
                ratio<-seq(input$simulate_ratios[2],input$simulate_ratios[1],(-1)*input$simulate_steps)
                ratio<-ratio/100
                thresholds_del<-data.frame(ratio=ratio,minSnps=NA,window=NA)
                thresholds_dup<-data.frame(ratio=ratio,minSnps=NA,window=NA)
                shinyjs::html("text", paste0("&nbsp&nbsp&nbsp&nbsp&nbsp&nbspA: Evaluating Tumor Ratios","<br>"), add = TRUE)  
                progress <- shiny::Progress$new()
                progress$set(message = "A: Evaluating Tumor Ratios", value = 0)
                for(j in 1:length(ratio)){
                    #message("Ratio: ",ratio[j])
                    progress$inc(input$simulate_steps/100)
                    thresholds_snps<-data.frame(snps=1:100,significant_del=0,significant_dup=0)
                    vaf_del<-(ratio[j]-1)/(-2+ratio[j])
                    vaf_g<-0.5  
                    vaf_dup<-1/(2+ratio[j])
                    simulations<-input$simulations_nr
                    
                    info_del<-data.frame(cov_tumor=rep(NA,100*simulations),vaf_tumor_del=rep(NA,100*simulations),
                                         vaf_tumor_lower_del=rep(NA,100*simulations),vaf_tumor_upper_del=rep(NA,100*simulations),
                                         cov_germline=rep(NA,100*simulations),vaf_germline=rep(NA,100*simulations),
                                         vaf_germline_lower=rep(NA,100*simulations),vaf_germline_upper=rep(NA,100*simulations),
                                         vaf_tumor_dup=rep(NA,100*simulations),
                                         vaf_tumor_lower_dup=rep(NA,100*simulations),vaf_tumor_upper_dup=rep(NA,100*simulations))
                    temp<-round(rlnorm(meanlog=mean(log(tumor[,8])),sdlog = sd(log(tumor[,8])),n=200*simulations))
                    info_del[,1]<-temp[temp>0][1:(100*simulations)]
                    temp<-round(rlnorm(meanlog=mean(log(germline[,8])),sdlog = sd(log(germline[,8])),n=200*simulations))
                    info_del[,5]<-temp[temp>0][1:(100*simulations)]
                    
                    info_del[,2]<-apply(info_del,MARGIN=1,
                                        FUN=function(x){sqrt((round(rnorm(mean=vaf_del*x[1],
                                                                          sd = 0.05*vaf_del*x[1],n=1)))**2)})/info_del[,1]
                    info_del[,3]<-apply(info_del,MARGIN=1,
                                        FUN=function(x){binom.test(round(x[1]*x[2]),
                                                                   x[1],x[2],conf.level = 0.95)$conf.int[1]})
                    info_del[,4]<-apply(info_del,MARGIN=1,
                                        FUN=function(x){binom.test(round(x[1]*x[2]),
                                                                   x[1],x[2],conf.level = 0.95)$conf.int[2]})
                    
                    info_del[,6]<-apply(info_del,MARGIN=1,
                                        FUN=function(x){sqrt((round(rnorm(mean=vaf_g*x[5],
                                                                          sd = 0.05*vaf_g*x[5],n=1)))**2)})/info_del[,5]
                    info_del[,7]<-apply(info_del,MARGIN=1,
                                        FUN=function(x){binom.test(round(x[5]*x[6]),
                                                                   x[5],x[6],conf.level = 0.95)$conf.int[1]})
                    info_del[,8]<-apply(info_del,MARGIN=1,
                                        FUN=function(x){binom.test(round(x[5]*x[6]),
                                                                   x[5],x[6],conf.level = 0.95)$conf.int[2]})
                    
                    info_del[,9]<-apply(info_del,MARGIN=1,
                                        FUN=function(x){sqrt((round(rnorm(mean=vaf_dup*x[1],
                                                                          sd = 0.05*vaf_dup*x[1],n=1)))**2)})/info_del[,1]
                    info_del[,10]<-apply(info_del,MARGIN=1,
                                         FUN=function(x){binom.test(round(x[1]*x[9]),
                                                                    x[1],x[9],conf.level = 0.95)$conf.int[1]})
                    info_del[,11]<-apply(info_del,MARGIN=1,
                                         FUN=function(x){binom.test(round(x[1]*x[9]),
                                                                    x[1],x[9],conf.level = 0.95)$conf.int[2]})   
                    
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
                    for(k in c(2,3,4)){
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
                    for(k in c(6,7,8)){
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
                    for(k in c(9,10,11)){
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
                    for(k in c(6,7,8)){
                        temp1<-info_del[,k]
                        temp1[info_del[,6]>0.5]<-NA
                        temp1.1<-(2*temp1-1)/(-1*temp1)
                        temp2<-info_del[,k]
                        temp2[info_del[,6]<=0.5]<-NA
                        temp2.1<-(2*temp2-1)/(1-temp2)
                        temp3<-apply(cbind(temp1.1,temp2.1),MARGIN=1,FUN=function(x){max(x,na.rm=T)})
                        info_del_r[,k+6]<-temp3
                    }
                    
                    temp<-apply(info_del_r,MARGIN = 1,
                                FUN=function(x){(abs(range(x[3],x[4])[2]-range(x[3],x[4])[1])*abs(range(x[7],x[8])[2]-range(x[7],x[8])[1]))**(-1)})
                    test_del<-cbind((info_del_r[,2]-info_del_r[,6]),temp)
                    temp<-apply(info_del_r,MARGIN = 1,
                                FUN=function(x){(abs(range(x[10],x[11])[2]-range(x[10],x[11])[1])*abs(range(x[13],x[14])[2]-range(x[13],x[14])[1]))**(-1)})
                    test_dup<-cbind((info_del_r[,9]-info_del_r[,12]),temp)

                    i<-2
                    del_enough<-F
                    dup_enough<-F
                    start<-seq(1,simulations*100,100)
                    while(i<=100&&(del_enough==F||dup_enough==F)){
                        for(k in 1:length(start)){
                            if(del_enough==F){
                                temp<-wtd.t.test(x=test_del[start[k]:(start[k]+i-1),1],
                                                 weight=test_del[start[k]:(start[k]+i-1),2],
                                                 alternative="greater")
                                if(!is.na(temp$coefficients[3])&&temp$coefficients[3]<0.025){
                                    thresholds_snps[i,2]<-sum(thresholds_snps[i,2],1)
                                } 
                            }
                            if(dup_enough==F){
                                temp<-wtd.t.test(x=test_dup[start[k]:(start[k]+i-1),1],
                                                 weight=test_dup[start[k]:(start[k]+i-1),2],
                                                 alternative="greater")
                                if(!is.na(temp$coefficients[3])&&temp$coefficients[3]<0.025){
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
                    thresholds_del[ratio[j]==thresholds_del[,1],2]<-thresholds_snps[min(which(thresholds_snps[,2]>=0.95),100,na.rm=T),1]
                    thresholds_dup[ratio[j]==thresholds_dup[,1],2]<-thresholds_snps[min(which(thresholds_snps[,3]>=0.95),100,na.rm=T),1]
                }
                progress$close()
                
                shinyjs::html("text", paste0("&nbsp&nbsp&nbsp&nbsp&nbsp&nbspB: Evaluating Window Sizes","<br>"), add = TRUE)  
                progress <- shiny::Progress$new()
                progress$set(message = "B: Evaluating Window Sizes", value = 0)
                chromosomes<-c(249250621,243199373,198022430,191154276,180915260,171115067,
                               159138663,146364022,141213431,135534747,135006516,133851895,
                               115169878,107349540,102531392,90354753,81195210,78077248,
                               59128983,63025520,48129895,51304566)
                
                snps_to_check<-unique(sort(c(thresholds_del[,2],thresholds_dup[,2])))
                window<-backup<-input$simulate_window_start
                results<-data.frame(Window=window,meanSnps=NA)
                for(investigate in snps_to_check){
                    progress$inc(1/(length(snps_to_check)+1))
                    window<-backup
                    while(window<=min(chromosomes)){
                        counter<-1
                        nr_snps<-data.frame(snps=seq(0,1000),nr=0)
                        for(chr in 1:22){
                            chr_data<-germline[germline[,1]==chr,c(1,2)]
                            start<-1
                            while((start+window)<=chromosomes[chr]){
                                temp<-chr_data[chr_data[,2]>=start&chr_data[,2]<(start+window),]
                                snps_in_window<-length(temp[,1])
                                first_snp<-temp[1,2]
                                next_snp<-chr_data[grep(paste("^",temp[length(temp[,1]),2],"$",sep=""),chr_data[,2])+1,2]
                                moved<-F
                                if(snps_in_window==0){
                                    nr_snps[1,2]<-nr_snps[1,2]+window
                                    start<-start+window
                                    moved<-T
                                }
                                if(moved==F&&is.na(next_snp)){
                                    next_snp<-chromosomes[chr]
                                }
                                if(moved==F&&!is.na(next_snp)&&length(next_snp)>1){
                                    next_snp<-next_snp[length(next_snp)]
                                }
                                if(moved==F&&(first_snp-start)>=(next_snp-(start+window))){
                                    nr_snps[nr_snps[,1]==snps_in_window,2]<-nr_snps[nr_snps[,1]==snps_in_window,2]+(next_snp-(start+window)+1)
                                    start<-next_snp-window+1
                                    moved<-T
                                }
                                if(moved==F&&(first_snp-start)<(next_snp-temp[length(temp[,1]),2])){
                                    nr_snps[nr_snps[,1]==snps_in_window,2]<-nr_snps[nr_snps[,1]==snps_in_window,2]+(first_snp-start)
                                    start<-first_snp+1
                                }
                            }
                        }
                        results[window==results[,1],2]<-sum(nr_snps[,1]*nr_snps[,2])/sum(nr_snps[,2])
                        if(max(results[,2],na.rm=T)>=investigate){
                            thresholds_del[thresholds_del[,2]==investigate,3]<-window
                            thresholds_dup[thresholds_dup[,2]==investigate,3]<-window
                            backup<-window
                            window<-max(chromosomes)
                        }
                        if(max(results[,2],na.rm=T)<investigate){
                            window<-window+input$simulate_window_steps
                            results<-rbind(results,data.frame(Window=window,meanSnps=NA))
                        }
                    }
                }
                progress$close()
                
                thresholds_del<-cbind(thresholds_del,windowTransformed=thresholds_del[,3]/max(thresholds_del[,3],na.rm=T))
                thresholds_dup<-cbind(thresholds_dup,windowTransformed=thresholds_dup[,3]/max(thresholds_dup[,3],na.rm=T))
                thresholds_del<-cbind(thresholds_del,Distance=sqrt(thresholds_del[,1]**2+thresholds_del[,4]**2))
                thresholds_dup<-cbind(thresholds_dup,Distance=sqrt(thresholds_dup[,1]**2+thresholds_dup[,4]**2))
                
                shinyjs::html("text", paste0("<br>","&nbsp&nbsp&nbspDetection thresholds for CNVs:","<br>"), add = TRUE)
                if(input$simulation_optimization=="Compromize (ratio and window size)"){
                    del_ratio<-thresholds_del[which.min(thresholds_del[,5]),1]
                    del_window<-thresholds_del[which.min(thresholds_del[,5]),3]
                    shinyjs::html("text", paste0("&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbspDeletions:",del_window,
                                                 " bp; ",del_ratio*100,"% Tumor cells","<br>"), add = TRUE)
                    
                    dup_ratio<-thresholds_dup[which.min(thresholds_dup[,5]),1]
                    dup_window<-thresholds_dup[which.min(thresholds_dup[,5]),3]
                    shinyjs::html("text", paste0("&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbspDuplications:",dup_window,
                                                 " bp; ",dup_ratio*100,"% Tumor cells","<br>"), add = TRUE)
                }
                if(input$simulation_optimization=="Force (ratio)"){
                    del_ratio<-thresholds_del[which.min(thresholds_del[!is.na(threshlds_del[,5]),1]),1]
                    del_window<-thresholds_del[which.min(thresholds_del[!is.na(threshlds_del[,5]),1]),3]
                    shinyjs::html("text", paste0("&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbspDeletions:",del_window,
                                                 " bp; ",del_ratio*100,"% Tumor cells","<br>"), add = TRUE)
                    
                    dup_ratio<-thresholds_dup[which.min(thresholds_dup[!is.na(threshlds_dup[,5]),1]),1]
                    dup_window<-thresholds_dup[which.min(thresholds_dup[!is.na(threshlds_dup[,5]),1]),3]
                    shinyjs::html("text", paste0("&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbspDuplications:",dup_window,
                                                 " bp; ",dup_ratio*100,"% Tumor cells","<br>"), add = TRUE)
                }
                if(input$simulation_optimization=="Force (window size)"){
                    del_ratio<-thresholds_del[which.min(thresholds_del[!is.na(threshlds_del[,5]),3]),1]
                    del_window<-thresholds_del[which.min(thresholds_del[!is.na(threshlds_del[,5]),3]),3]
                    shinyjs::html("text", paste0("&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbspDeletions:",del_window,
                                                 " bp; ",del_ratio*100,"% Tumor cells","<br>"), add = TRUE)
                    
                    dup_ratio<-thresholds_dup[which.min(thresholds_dup[!is.na(threshlds_dup[,5]),3]),1]
                    dup_window<-thresholds_dup[which.min(thresholds_dup[!is.na(threshlds_dup[,5]),3]),3]
                    shinyjs::html("text", paste0("&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbspDuplications:",dup_window,
                                                 " bp; ",dup_ratio*100,"% Tumor cells","<br>"), add = TRUE)
                }
                
                if(n==1){
                    detectionThresholds<-data.frame(Sample=c(samples_t[n,1]),Window_del=c(del_window),
                                                    Ratio_del=c(del_ratio),Window_dup=c(dup_window),
                                                    Ratio_dup=c(dup_ratio))
                }
                if(n>1){
                    detectionThresholds<-rbind(detectionThresholds,data.frame(Sample=c(samples_t[n,1]),
                                                                              Window_del=c(del_window),
                                                                              Ratio_del=c(del_ratio),
                                                                              Window_dup=c(dup_window),
                                                                              Ratio_dup=c(dup_ratio)))
                }
                output$table_dt <- renderDataTable(datatable(detectionThresholds))
                write.table(detectionThresholds,paste0(input$output_folder,"DetectionThresholds.txt"),
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
            counter<-1
            progress <- shiny::Progress$new()
            progress$set(message = "Detecting deletions", value = 0)
            for(chr in 1:22){
                progress$inc(1/23)
                tumor_1<-tumor[tumor[,1]==chr,]
                del_1<-del[tumor[,1]==chr,]
                del_g1<-del_g[tumor[,1]==chr,]
                for(i in 0:(length(tumor_1[,1])-1)){
                    if(i==0){
                        start<-1
                    }
                    if(i>0){
                        start<-tumor_1[i,2]+1
                    }
                    
                    del_oi<-del_1[tumor_1[,2]>=start&tumor_1[,2]<=(start+window),]
                    del_goi<-del_g1[tumor_1[,2]>=start&tumor_1[,2]<=(start+window),]
                    if(length(del_oi[,1])>0&&length(del_goi[,1])>0&&
                       sum(!is.na(del_oi[,2]))>1&&sum(!is.na(del_goi[,2]))>1&&
                       sum(!is.na(del_oi[,2]-del_goi[,2]))>1&&
                       sum(!is.na(unique(round(del_oi[,2]-del_goi[,2],digits = 7))))>1){
                        temp<-wtd.t.test(x=del_oi[,2]-del_goi[,2],
                                         weight=(abs(del_oi[,1]-del_oi[,3])*abs(del_goi[,1]-del_goi[,3]))**(-1),
                                         alternative = "greater")
                        temp2<-wtd.t.test(x=del_oi[,2]-del_goi[,2],
                                          weight=(abs(del_oi[,1]-del_oi[,3])*abs(del_goi[,1]-del_goi[,3]))**(-1),
                                          alternative = "less")
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
                        counter<-counter+1
                    }
                }
            } 
            progress$close()
            
            #Calling duplications:
            window<-detectionThresholds[n,4]
            progress <- shiny::Progress$new()
            progress$set(message = "Detecting duplications", value = 0)
            for(chr in 1:22){
                progress$inc(1/23)
                tumor_1<-tumor[tumor[,1]==chr,]
                dup_1<-dup[tumor[,1]==chr,]
                dup_g1<-dup_g[tumor[,1]==chr,]
                for(i in 0:(length(tumor_1[,1])-1)){
                    if(i==0){
                        start<-1
                    }
                    if(i>0){
                        start<-tumor_1[i,2]+1
                    }
                    dup_oi<-dup_1[tumor_1[,2]>=start&tumor_1[,2]<=(start+window),]
                    dup_goi<-dup_g1[tumor_1[,2]>=start&tumor_1[,2]<=(start+window),]
                    
                    if(length(dup_oi[,1])>0&&length(dup_goi[,1])>0&&
                       sum(!is.na(dup_oi[,2]))>1&&sum(!is.na(dup_goi[,2]))>1&&
                       sum(!is.na(dup_oi[,2]-dup_goi[,2]))>1&&
                       sum(!is.na(unique(round(dup_oi[,2]-dup_goi[,2],digits = 7))))>1){
                        temp<-wtd.t.test(x=dup_oi[,2]-dup_goi[,2],
                                         weight=(abs(dup_oi[,1]-dup_oi[,3])*abs(dup_goi[,1]-dup_goi[,3]))**(-1),
                                         alternative = "greater")
                        temp2<-wtd.t.test(x=dup_oi[,2]-dup_goi[,2],
                                          weight=(abs(dup_oi[,1]-dup_oi[,3])*abs(dup_goi[,1]-dup_goi[,3]))**(-1),
                                          alternative = "less")
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
                        counter<-counter+1
                    }
                }
            }
            progress$close()
            
            results<-results[order(results[,4],results[,1],results[,2]),]
            results20<-results[results[,8]>=20&((results[,4]=="del"&results[,6]>=detectionThresholds[n,3])|(results[,4]=="dup"&results[,6]>=detectionThresholds[n,5])),]
            if(sum(input$output_files=="Raw CNV calls")>0){
                write.table(results20,paste(input$output_folder,samples_t[n,1],".CNVs_raw.txt",sep=""),row.names=F,sep="\t",quote=F)
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
                
                png(paste0(input$output_folder,samples_t[n,1],"_raw_all.png"),width=1800,height=800)
                plot(NULL,xlim=c(0,2881033286),ylim=c(-1,1),xaxt="n",yaxt="n",
                     xlab="Choromosome",ylab="",main=paste("Sample ",samples_t[n,1],sep=""))
                colfunc <- colorRampPalette(c("blue","red"))
                colpalette<-colfunc((max(x.value_del[,3])-min(x.value_del[,3]))*100+1)
                for(i in 1:length(x.value_del[,1])){
                    points(x.value_del[i,c(1,2)],c(x.value_del[i,4],x.value_del[i,4]),type="l",lwd=5,
                           col=colpalette[min(round(100*x.value_del[i,3])+1,length(colpalette))])
                }
                colpalette<-colfunc((max(x.value_dup[,3])-min(x.value_dup[,3]))*100+1)
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
                    list(src=paste0(input$output_folder,samples_t[n,1],"_raw_all.png"),
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
                    
                    png(paste0(input$output_folder,samples_t[n,1],"_raw_sig.png"),width=1800,height=800)
                    plot(NULL,xlim=c(0,2881033286),ylim=c(-1,1),xaxt="n",yaxt="n",
                         xlab="Choromosome",ylab="",main=paste("Sample ",samples_t[n,1],sep=""))
                    colfunc <- colorRampPalette(c("blue","red"))
                    if(no_dels==F){
                        colpalette<-colfunc((max(x.value_del[,3])-min(x.value_del[,3]))*100+1)
                        for(i in 1:length(x.value_del[,1])){
                            points(x.value_del[i,c(1,2)],c(x.value_del[i,4],x.value_del[i,4]),type="l",lwd=5,
                                   col=colpalette[min(round(100*x.value_del[i,3])+1,length(colpalette))])
                        }
                    }
                    if(no_dups==F){
                        colpalette<-colfunc((max(x.value_dup[,3])-min(x.value_dup[,3]))*100+1)
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
                        list(src=paste0(input$output_folder,samples_t[n,1],"_raw_sig.png"),
                             height=400,
                             width=900)},
                        deleteFile = FALSE
                    )  
                }
                
                results20<-cbind(results20,windows=1)
                window2<-input$maxDist
                results_merged<-data.frame(chr=NA,start=NA,end=NA,variant=NA,p.value_min=NA,
                                           p.value_max=NA,cells_min=NA,cells_max=NA,
                                           cells_avg=NA,snps=NA,qual=NA,sd_low=NA,
                                           sd_high=NA,windows=NA)
                meanList<-list()
                counter<-1
                
                if(length(results20[,1])>0&&!is.na(results20[1,1])){
                    shinyjs::html("text", paste0("&nbsp&nbsp&nbspMerging raw CNVs","<br>"), add = TRUE)
                    progress_sample$inc(1/4)
                    for(i in 1:(length(results20[,1]))){
                        if(i==1){
                            results_merged[counter,1]<-results20[i,1]
                            results_merged[counter,2]<-results20[i,2]
                            results_merged[counter,3]<-results20[i,3]
                            results_merged[counter,4]<-results20[i,4]
                            results_merged[counter,5]<-results20[i,5]
                            results_merged[counter,6]<-results20[i,5]
                            results_merged[counter,7]<-results20[i,6]
                            results_merged[counter,8]<-results20[i,6]
                            results_merged[counter,10]<-results20[i,7]
                            results_merged[counter,11]<-results20[i,8]
                            temp_avg<-results20[i,6]
                            results_merged[counter,12]<-temp_sd_low<-results20[i,6]-2*results20[i,9]
                            results_merged[counter,13]<-temp_sd_high<-results20[i,6]+2*results20[i,9]
                            results_merged[counter,14]<-results20[i,10]
                        }
                        if(i>1){
                            if(results_merged[counter,3]>=results20[i,2]&&
                               results_merged[counter,4]==results20[i,4]&&
                               results_merged[counter,1]==results20[i,1]&&
                               results20[i,6]>=temp_sd_low&&results20[i,6]<=temp_sd_high){
                                #overlapping; same variant -> merge
                                results_merged[counter,3]<-results20[i,3]
                                if(results20[i,5]<results_merged[counter,5]){
                                    results_merged[counter,5]<-results20[i,5]
                                }
                                if(results20[i,5]>results_merged[counter,6]){
                                    results_merged[counter,6]<-results20[i,5]
                                }
                                if(results20[i,6]<results_merged[counter,7]){
                                    results_merged[counter,7]<-results20[i,6]
                                }
                                if(results20[i,6]>results_merged[counter,8]){
                                    results_merged[counter,8]<-results20[i,6]
                                }
                                if((results20[i,6]-results20[i,9])<temp_sd_low){
                                    results_merged[counter,12]<-temp_sd_low<-results20[i,6]-2*results20[i,9]
                                }
                                if((results20[i,6]+results20[i,9])>temp_sd_high){
                                    results_merged[counter,13]<-temp_sd_high<-results20[i,6]+2*results20[i,9]
                                }
                                temp_sd_low<-temp_sd_low+0.05*temp_sd_low
                                temp_sd_high<-temp_sd_high-0.05*temp_sd_high
                                results_merged[counter,10]<-sum(results_merged[counter,10],results20[i,7],na.rm=T)   
                                results_merged[counter,11]<-sum(results_merged[counter,11],results20[i,8],na.rm=T)
                                results_merged[counter,14]<-sum(results_merged[counter,14],results20[i,10])
                                temp_avg<-c(temp_avg,results20[i,6])
                            }
                            if(results_merged[counter,3]<results20[i,2]&&
                               results_merged[counter,4]==results20[i,4]&&
                               results_merged[counter,1]==results20[i,1]&&
                               results20[i,6]>=temp_sd_low&&results20[i,6]<=temp_sd_high){
                                #not overlapping, but same variant -> maybe merge
                                start_line<-(which(results20[i-1,3]==results[,3]&results20[i-1,1]==results[,1]&results20[i-1,4]==results[,4])+1)[1]
                                end_line<-(which(results20[i,2]==results[,2]&results20[i,1]==results[,1]&results20[i,4]==results[,4])-1)[1]
                                if(end_line<start_line){
                                    #merging unproblematic; no regions in between
                                    results_merged[counter,3]<-results20[i,3]
                                    if(results20[i,5]<results_merged[counter,5]){
                                        results_merged[counter,5]<-results20[i,5]
                                    }
                                    if(results20[i,5]>results_merged[counter,6]){
                                        results_merged[counter,6]<-results20[i,5]
                                    }
                                    if(results20[i,6]<results_merged[counter,7]){
                                        results_merged[counter,7]<-results20[i,6]
                                    }
                                    if(results20[i,6]>results_merged[counter,8]){
                                        results_merged[counter,8]<-results20[i,6]
                                    }
                                    if((results20[i,6]-results20[i,9])<temp_sd_low){
                                        results_merged[counter,12]<-temp_sd_low<-results20[i,6]-2*results20[i,9]
                                    }
                                    if((results20[i,6]+results20[i,9])>temp_sd_high){
                                        results_merged[counter,13]<-temp_sd_high<-results20[i,6]+2*results20[i,9]
                                    }
                                    temp_sd_low<-temp_sd_low+0.05*temp_sd_low
                                    temp_sd_high<-temp_sd_high-0.05*temp_sd_high
                                    results_merged[counter,10]<-sum(results_merged[counter,10],results20[i,7],na.rm=T)   
                                    results_merged[counter,11]<-sum(results_merged[counter,11],results20[i,8],na.rm=T)
                                    results_merged[counter,14]<-sum(results_merged[counter,14],results20[i,10])
                                    temp_avg<-c(temp_avg,results20[i,6])
                                }
                                if(end_line>=start_line){
                                    #evaluate regions in between, then: decide if merge
                                    problematic<-results[start_line:end_line,]
                                    if(sum(problematic[,8]<=-20)==0&&(results20[i,2]-results20[i-1,3])<window2){
                                        #no region significant for "less" and not too far away -> merge
                                        results_merged[counter,3]<-results20[i,3]
                                        if(results20[i,5]<results_merged[counter,5]){
                                            results_merged[counter,5]<-results20[i,5]
                                        }
                                        if(results20[i,5]>results_merged[counter,6]){
                                            results_merged[counter,6]<-results20[i,5]
                                        }
                                        if(results20[i,6]<results_merged[counter,7]){
                                            results_merged[counter,7]<-results20[i,6]
                                        }
                                        if(results20[i,6]>results_merged[counter,8]){
                                            results_merged[counter,8]<-results20[i,6]
                                        }
                                        if((results20[i,6]-results20[i,9])<temp_sd_low){
                                            results_merged[counter,12]<-temp_sd_low<-results20[i,6]-2*results20[i,9]
                                        }
                                        if((results20[i,6]+results20[i,9])>temp_sd_high){
                                            results_merged[counter,13]<-temp_sd_high<-results20[i,6]+2*results20[i,9]
                                        }
                                        temp_sd_low<-temp_sd_low+0.05*temp_sd_low
                                        temp_sd_high<-temp_sd_high-0.05*temp_sd_high
                                        results_merged[counter,10]<-sum(results_merged[counter,10],results20[i,7],na.rm=T)   
                                        results_merged[counter,11]<-sum(results_merged[counter,11],results20[i,8],na.rm=T)
                                        results_merged[counter,14]<-sum(results_merged[counter,14],results20[i,10])
                                        temp_avg<-c(temp_avg,results20[i,6])
                                    }
                                    if(sum(problematic[,8]<=-20)>0||(results20[i,2]-results20[i-1,3])>=window2){
                                        #don't merge
                                        results_merged[counter,9]<-mean(temp_avg)
                                        meanList[[counter]]<-temp_avg
                                        counter<-counter+1
                                        results_merged[counter,1]<-results20[i,1]
                                        results_merged[counter,2]<-results20[i,2]
                                        results_merged[counter,3]<-results20[i,3]
                                        results_merged[counter,4]<-results20[i,4]
                                        results_merged[counter,5]<-results20[i,5]
                                        results_merged[counter,6]<-results20[i,5]
                                        results_merged[counter,7]<-results20[i,6]
                                        results_merged[counter,8]<-results20[i,6]
                                        results_merged[counter,10]<-results20[i,7]
                                        results_merged[counter,11]<-results20[i,8]
                                        temp_avg<-results20[i,6]
                                        results_merged[counter,12]<-temp_sd_low<-results20[i,6]-2*results20[i,9]
                                        results_merged[counter,13]<-temp_sd_high<-results20[i,6]+2*results20[i,9]
                                        results_merged[counter,14]<-results20[i,10]
                                    }
                                }
                            }
                            if(results_merged[counter,4]!=results20[i,4]||
                               results_merged[counter,1]!=results20[i,1]||
                               results20[i,6]<temp_sd_low||results20[i,6]>temp_sd_high){
                                #don't merge
                                results_merged[counter,9]<-mean(temp_avg)
                                meanList[[counter]]<-temp_avg
                                counter<-counter+1
                                results_merged[counter,1]<-results20[i,1]
                                results_merged[counter,2]<-results20[i,2]
                                results_merged[counter,3]<-results20[i,3]
                                results_merged[counter,4]<-results20[i,4]
                                results_merged[counter,5]<-results20[i,5]
                                results_merged[counter,6]<-results20[i,5]
                                results_merged[counter,7]<-results20[i,6]
                                results_merged[counter,8]<-results20[i,6]
                                results_merged[counter,10]<-results20[i,7]
                                results_merged[counter,11]<-results20[i,8]
                                temp_avg<-results20[i,6]
                                results_merged[counter,12]<-temp_sd_low<-results20[i,6]-2*results20[i,9]
                                results_merged[counter,13]<-temp_sd_high<-results20[i,6]+2*results20[i,9]
                                results_merged[counter,14]<-results20[i,10]
                            }
                        }
                    }
                    results_merged[counter,9]<-mean(temp_avg)
                    meanList[[counter]]<-temp_avg
                    
                    for(i in 1:length(results_merged[,1])){
                        if(!is.na(results_merged[i,1])){
                            for(j in i+1:length(results_merged[,1])){
                                if(!is.na(results_merged[j,1])&&
                                   results_merged[i,1]==results_merged[j,1]&&
                                   results_merged[i,4]==results_merged[j,4]&&
                                   results_merged[i,12]<=results_merged[j,9]&&
                                   results_merged[i,13]>=results_merged[j,9]&&
                                   (results_merged[i,3]+window2)>=results_merged[j,2]){
                                    results_merged[i,3]<-max(results_merged[i,3],results_merged[j,3])
                                    if(results_merged[j,5]<results_merged[i,5]){
                                        results_merged[i,5]<-results_merged[j,5]
                                    }
                                    if(results_merged[j,6]>results_merged[i,6]){
                                        results_merged[i,6]<-results_merged[j,6]
                                    }
                                    if(results_merged[j,7]<results_merged[i,7]){
                                        results_merged[i,7]<-results_merged[j,7]
                                    }
                                    if(results_merged[j,8]>results_merged[i,8]){
                                        results_merged[i,8]<-results_merged[j,8]
                                    }
                                    results_merged[i,9]<-mean(c(meanList[[i]],meanList[[j]]))
                                    meanList[[i]]<-c(meanList[[i]],meanList[[j]])
                                    results_merged[i,10]<-sum(results_merged[i,10],results_merged[j,10],na.rm=T)
                                    results_merged[i,11]<-sum(results_merged[i,11],results_merged[j,11],na.rm=T)
                                    results_merged[i,14]<-sum(results_merged[i,14],results_merged[j,14],na.rm=T)
                                    if(results_merged[j,12]<results_merged[i,12]){
                                        results_merged[i,12]<-results_merged[j,12]
                                    }
                                    if(results_merged[j,13]>results_merged[i,13]){
                                        results_merged[i,13]<-results_merged[j,13]
                                    }
                                    temp_length<-0.025*(results_merged[i,13]-results_merged[i,12])
                                    results_merged[i,12]<-results_merged[i,12]+temp_length
                                    results_merged[i,13]<-results_merged[i,13]-temp_length
                                    results_merged[j,1]<-NA
                                }
                            }     
                        }
                    }
                    results_merged<-results_merged[!is.na(results_merged[,1]),]
                    results_merged<-results_merged[results_merged[,9]<1.5,]
                    results_merged<-cbind(results_merged,logQual=log(results_merged[,11]))
                    results_merged[,5]<-round(results_merged[,5],digits = 4)
                    results_merged[,6]<-round(results_merged[,6],digits = 4)
                    results_merged[,7]<-round(results_merged[,7],digits = 2)
                    results_merged[,8]<-round(results_merged[,8],digits = 2)
                    results_merged[,9]<-round(results_merged[,9],digits = 2)
                    results_merged[,15]<-round(results_merged[,15],digits = 2)
                    results_merged<-cbind(Sample=samples_t[n,1],results_merged,stringsAsFactors=F)
                    results_merged_report<-results_merged[,c(1:11,15,16)]
                    results_merged_report[,13]<-as.character(results_merged_report[,13])
                    if(sum(input$output_files=="Merged CNV calls")>0){
                        write.table(results_merged,paste0(input$output_folder,samples_t[n,1],".CNVs_merged.txt"),row.names=F,sep="\t",quote=F)   
                    }
                    if(sum(input$output_plots=="Merged CNV calls")>0){
                        output$text_plot3<-renderText({"Merged CNV calls"})
                        #plot
                        no_dels<-T
                        x.value_del<-data.frame(start=NA,end=NA,quality=NA,cells=NA)
                        results_merged_del<-results_merged[results_merged[,5]=="del",c(2:16)]
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
                                x.value_del[i,3]<-log(results_merged_del[i,11])
                                x.value_del[i,4]<-results_merged_del[i,9]*(-1)
                            }
                        }
                        if(length(x.value_del[,1])>1||!is.na(x.value_del[,1])){
                            x.value_del[x.value_del[,3]==Inf|x.value_del[,3]==-Inf,3]<-NA
                            x.value_del[is.na(x.value_del[,3]),3]<-min(max(x.value_del[,3],1,na.rm=T),100,na.rm=T)
                            no_dels<-F
                        }
                        no_dups<-T
                        x.value_dup<-data.frame(start=NA,end=NA,quality=NA,cells=NA)
                        results_merged_dup<-results_merged[results_merged[,5]=="dup",c(2:16)]
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
                                x.value_dup[i,3]<-log(results_merged_dup[i,11])
                                x.value_dup[i,4]<-results_merged_dup[i,9]
                            }
                        }
                        if(length(x.value_dup[,1])>1||!is.na(x.value_dup[,1])){
                            x.value_dup[x.value_dup[,3]==Inf|x.value_dup[,3]==-Inf,3]<-NA
                            x.value_dup[is.na(x.value_dup[,3]),3]<-min(max(x.value_dup[,3],1,na.rm=T),100,na.rm=T) 
                            no_dups<-F
                        }
                        png(paste0(input$output_folder,samples_t[n,1],"_merged.png"),width=1800,height=800)
                        plot(NULL,xlim=c(0,2881033286),ylim=c(-1,1),xaxt="n",yaxt="n",
                             xlab="Choromosome",ylab="",main=paste("Sample ",samples_t[n,1],sep=""))
                        colfunc <- colorRampPalette(c("blue","red"))
                        if(no_dels==F){
                            colpalette<-colfunc((max(x.value_del[,3])-min(x.value_del[,3]))*100+1)
                            for(i in 1:length(x.value_del[,1])){
                                points(x.value_del[i,c(1,2)],c(x.value_del[i,4],x.value_del[i,4]),type="l",lwd=5,
                                       col=colpalette[min(round(100*x.value_del[i,3])+1,length(colpalette))])
                            }            
                        }
                        if(no_dups==F){
                            colpalette<-colfunc((max(x.value_dup[,3])-min(x.value_dup[,3]))*100+1)
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
                            list(src=paste0(input$output_folder,samples_t[n,1],"_merged.png"),
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
                        results_final<-results_merged[results_merged[,16]>=input$quality_filter,]
                        results_final_report<-results_final[,c(1:11,15,16)]
                        results_final_report[,13]<-as.character(results_final_report[,13])
                        output$table_cnvs <- renderDataTable(datatable(results_final_report))
                        if(sum(input$output_files=="Filtered CNV calls")>0){
                            write.table(results_final_report,paste0(input$output_folder,samples_t[n,1],".CNVs_filtered.txt"),row.names=F,sep="\t",quote=F)   
                        }
                        
                        if(sum(input$output_plots=="Filtered CNV calls")>0){
                            output$text_plot4<-renderText({"Filtered CNV calls"})
                            #plot
                            no_dels<-T
                            x.value_del<-data.frame(start=NA,end=NA,quality=NA,cells=NA)
                            results_final_del<-results_final[results_final[,5]=="del",c(2:16)]
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
                                    x.value_del[i,3]<-log(results_final_del[i,11])
                                    x.value_del[i,4]<-results_final_del[i,9]*(-1)
                                }
                            }
                            if(length(x.value_del[,1])>1||!is.na(x.value_del[,1])){
                                x.value_del[x.value_del[,3]==Inf|x.value_del[,3]==-Inf,3]<-NA
                                x.value_del[is.na(x.value_del[,3]),3]<-min(max(x.value_del[,3],1,na.rm=T),100,na.rm=T)
                                no_dels<-F
                            }
                            no_dups<-T
                            x.value_dup<-data.frame(start=NA,end=NA,quality=NA,cells=NA)
                            results_final_dup<-results_final[results_final[,5]=="dup",c(2:16)]
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
                                    x.value_dup[i,3]<-log(results_final_dup[i,11])
                                    x.value_dup[i,4]<-results_final_dup[i,9]
                                }
                            }
                            if(length(x.value_dup[,1])>1||!is.na(x.value_dup[,1])){
                                x.value_dup[x.value_dup[,3]==Inf|x.value_dup[,3]==-Inf,3]<-NA
                                x.value_dup[is.na(x.value_dup[,3]),3]<-min(max(x.value_dup[,3],1,na.rm=T),100,na.rm=T)
                                no_dups<-F
                            }
                            
                            png(paste0(input$output_folder,samples_t[n,1],"_filtered.png"),width=1800,height=800)
                            plot(NULL,xlim=c(0,2881033286),ylim=c(-1,1),xaxt="n",yaxt="n",
                                 xlab="Choromosome",ylab="",main=paste("Sample ",samples_t[n,1],sep=""))
                            colfunc <- colorRampPalette(c("blue","red"))
                            if(no_dels==F){
                                colpalette<-colfunc((max(x.value_del[,3])-min(x.value_del[,3]))*100+1)
                                for(i in 1:length(x.value_del[,1])){
                                    points(x.value_del[i,c(1,2)],c(x.value_del[i,4],x.value_del[i,4]),type="l",lwd=5,
                                           col=colpalette[min(round(100*x.value_del[i,3])+1,length(colpalette))])
                                }            
                            }
                            if(no_dups==F){
                                colpalette<-colfunc((max(x.value_dup[,3])-min(x.value_dup[,3]))*100+1)
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
                                list(src=paste0(input$output_folder,samples_t[n,1],"_filtered.png"),
                                     height=400,
                                     width=900)},
                                deleteFile = FALSE
                            )  
                        }
                        
                    }
                }
            }
           
            
            if(length(results20[,1])==0){
                shinyjs::html("text", paste0("<br>","&nbsp&nbsp&nbspNo CNVs detected","<br>"), add = TRUE)
            }
            progress_sample$close()
        }
        updateRadioButtons(session,"select_samples",choices=samples_t[,1],
                           selected = samples_t[length(samples_t[,1]),1],inline=T)
        
        write.table(detectionThresholds,paste0(input$output_folder,"DetectionThresholds.txt"),
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
               file.exists(paste0(input$output_folder,input$select_samples,"_raw_all.png"))){
                output$text_plot1<-renderText({"Raw CNV calls (all)"})
                output$plot1 <- renderImage({
                    list(src=paste0(input$output_folder,input$select_samples,"_raw_all.png"),
                         height=400,
                         width=900)},
                    deleteFile = FALSE
                )  
            }
            if(sum(input$output_plots2=="Raw CNV calls (sig)")>0&&
               file.exists(paste0(input$output_folder,input$select_samples,"_raw_sig.png"))){
                if(sum(input$output_plots2=="Raw CNV calls (all)")>0){
                    output$text_plot2<-renderText({"Raw CNV calls (sig)"})
                    output$plot2 <- renderImage({
                        list(src=paste0(input$output_folder,input$select_samples,"_raw_sig.png"),
                             height=400,
                             width=900)},
                        deleteFile = FALSE
                    )  
                }
                if(sum(input$output_plots2=="Raw CNV calls (all)")==0){
                    output$text_plot1<-renderText({"Raw CNV calls (sig)"})
                    output$plot1 <- renderImage({
                        list(src=paste0(input$output_folder,input$select_samples,"_raw_sig.png"),
                             height=400,
                             width=900)},
                        deleteFile = FALSE
                    )  
                }

            }
            if(sum(input$output_plots2=="Merged CNV calls")>0&&
               file.exists(paste0(input$output_folder,input$select_samples,"_merged.png"))){
                if(sum(input$output_plots2=="Raw CNV calls (all)")>0){
                    if(sum(input$output_plots2=="Raw CNV calls (sig)")>0){
                        output$text_plot3<-renderText({"Merged CNV calls"})
                        output$plot3 <- renderImage({
                            list(src=paste0(input$output_folder,input$select_samples,"_merged.png"),
                                 height=400,
                                 width=900)},
                            deleteFile = FALSE
                        ) 
                    }
                    if(sum(input$output_plots2=="Raw CNV calls (sig)")==0){
                        output$text_plot2<-renderText({"Merged CNV calls"})
                        output$plot2 <- renderImage({
                            list(src=paste0(input$output_folder,input$select_samples,"_merged.png"),
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
                            list(src=paste0(input$output_folder,input$select_samples,"_merged.png"),
                                 height=400,
                                 width=900)},
                            deleteFile = FALSE
                        ) 
                    }
                    if(sum(input$output_plots2=="Raw CNV calls (sig)")==0){
                        output$text_plot1<-renderText({"Merged CNV calls"})
                        output$plot1 <- renderImage({
                            list(src=paste0(input$output_folder,input$select_samples,"_merged.png"),
                                 height=400,
                                 width=900)},
                            deleteFile = FALSE
                        ) 
                    }
                }
 
            }
            if(sum(input$output_plots2=="Filtered CNV calls")>0&&
               file.exists(paste0(input$output_folder,input$select_samples,"_filtered.png"))){
                if(sum(input$output_plots2=="Raw CNV calls (all)")>0){
                    if(sum(input$output_plots2=="Raw CNV calls (sig)")>0){
                        if(sum(input$output_plots2=="Merged CNV calls")>0){
                            output$text_plot4<-renderText({"Filtered CNV calls"})
                            output$plot4 <- renderImage({
                                list(src=paste0(input$output_folder,input$select_samples,"_filtered.png"),
                                     height=400,
                                     width=900)},
                                deleteFile = FALSE
                            ) 
                        }
                        if(sum(input$output_plots2=="Merged CNV calls")==0){
                            output$text_plot3<-renderText({"Filtered CNV calls"})
                            output$plot3 <- renderImage({
                                list(src=paste0(input$output_folder,input$select_samples,"_filtered.png"),
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
                                list(src=paste0(input$output_folder,input$select_samples,"_filtered.png"),
                                     height=400,
                                     width=900)},
                                deleteFile = FALSE
                            ) 
                        }
                        if(sum(input$output_plots2=="Merged CNV calls")==0){
                            output$text_plot2<-renderText({"Filtered CNV calls"})
                            output$plot2 <- renderImage({
                                list(src=paste0(input$output_folder,input$select_samples,"_filtered.png"),
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
                                list(src=paste0(input$output_folder,input$select_samples,"_filtered.png"),
                                     height=400,
                                     width=900)},
                                deleteFile = FALSE
                            ) 
                        }
                        if(sum(input$output_plots2=="Merged CNV calls")==0){
                            output$text_plot2<-renderText({"Filtered CNV calls"})
                            output$plot2 <- renderImage({
                                list(src=paste0(input$output_folder,input$select_samples,"_filtered.png"),
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
                                list(src=paste0(input$output_folder,input$select_samples,"_filtered.png"),
                                     height=400,
                                     width=900)},
                                deleteFile = FALSE
                            ) 
                        }
                        if(sum(input$output_plots2=="Merged CNV calls")==0){
                            output$text_plot1<-renderText({"Filtered CNV calls"})
                            output$plot1 <- renderImage({
                                list(src=paste0(input$output_folder,input$select_samples,"_filtered.png"),
                                     height=400,
                                     width=900)},
                                deleteFile = FALSE
                            ) 
                        }
                    }
                }
 
            }
            
            if(input$output_files2=="Merged CNV calls"&&
               file.exists(paste0(input$output_folder,input$select_samples,".CNVs_merged.txt"))){
                calls<-read.table(paste0(input$output_folder,input$select_samples,".CNVs_merged.txt"),
                                  header=T,quote = "",sep="\t",stringsAsFactors = F)
                calls[,13]<-as.character(calls[,13])
                output$table_cnvs <- renderDataTable(datatable(calls))
            }
            if(input$output_files2=="Filtered CNV calls"&&
               file.exists(paste0(input$output_folder,input$select_samples,".CNVs_filtered.txt"))){
                calls<-read.table(paste0(input$output_folder,input$select_samples,".CNVs_filtered.txt"),
                                  header=T,quote = "",sep="\t",stringsAsFactors = F)
                calls[,13]<-as.character(calls[,13])
                output$table_cnvs <- renderDataTable(datatable(calls))
            }
        }

    })
}
)
