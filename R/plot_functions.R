#' sr_plot function
#'
#' This function generates a stock-recruitment (S-R) model for rstan based on model inputs.
#' @param type Specify whether to generate a 'static' S-R model, where parameters are time-invariant, 
#' a time-varying 'tv' model, or a regime shift model 'regime'
#' @param df stock-recruitment dataset. With columns R (recruits) and S (spawners) and by (brood year)
#' @param form Either fit with Stan ('stan') or TMB ('tmb')
#' @param mod a fitted Stan or TMB model
#' @param pdf whether to create a pdf or not
#' @return returns the compiled rstan code for a given S-R model
#' @importFrom rstan stan_model
#' @export
#' @examples
#' sr_plot(type='static',df=df,form='stan',df=df,mod=f1,pdf=FALSE)

sr_plot=function(df,mod,title,make.pdf=FALSE,path,type=c('static','rw','hmm'),par=c('a','b','both'),form=c('stan','tmb'),ac=FALSE,sr_only=FALSE){
  if(type=='static'){ #static====
    x_new=seq(0,max(df$S),length.out=200)

    if(form=='stan'){
      post=rstan::extract(mod)
      pred_df=data.frame(pred=exp(median(post$log_a)-median(post$b)*x_new)*x_new,x_new=x_new)
      }
    if(form=='tmb'){
      pred_df=data.frame(pred=exp(mod$alpha-mod$beta*x_new)*x_new,x_new=x_new)
    }
    plot=ggplot2::ggplot(df, aes(S, R)) +
      geom_line(data=pred_df,aes(x=x_new,y=pred),linewidth=1.3)+
      geom_point(aes(colour = by),size=2.5) +
      scale_colour_viridis_c(name='Year')+
      ggtitle(title)+
      xlab("Spawners") + 
      ylab("Recruits")+
      xlim(0, max(df$S))+
      ylim(0, max(df$R))+
      theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
            strip.text = element_text(face="bold", size=12),
            axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
  }
  if(type=='rw'){
    x_new=seq(min(df$S),max(df$S),length.out=200)
    by_q=round(quantile(df$by,seq(0,1,by=0.1)))
      if(par=='a'){ #rw alpha=====
        if(form=='stan'){
          post=rstan::extract(mod)
          pred_df=data.frame(x_new)
        for(n in 1:length(by_q)){
          pred_df[,1+n]=exp(median(post$log_a[,match(by_q[n],df$by)])-median(post$b)*x_new)*x_new
        }
          alpha_df=data.frame(by=seq(min(df$by),max(df$by)),med=apply(post$log_a,2,median),l90=apply(post$log_a,2,quantile,0.1),u90=apply(post$log_a,2,quantile,0.9))
          plot2=ggplot2::ggplot(alpha_df, aes(by,med)) +
            geom_line(aes(x=by,y=med),linewidth=1.3)+
            geom_point(aes(colour = by),size=4) +
            scale_colour_viridis_c(name='Year')+
            geom_ribbon(aes(ymin =l90, ymax =u90), alpha = 0.2)+
            ggtitle(title)+
            xlab("Year") + 
            ylab("log(Alpha)")+
            theme_classic(14)+
            theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                  strip.text = element_text(face="bold", size=12),
                  axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
          
        }
        
        if(form=='tmb'){
          pred_df=data.frame(x_new)
          for(n in 1:length(by_q)){
            pred_df[,1+n]=exp(mod$alpha[match(by_q[n],df$by)]-mod$beta*x_new)*x_new
          }
          alpha_df=data.frame(by=df$by,med=mod$alpha)
          
          plot2=ggplot2::ggplot(alpha_df, aes(by,med)) +
            geom_hline(yintercept=0,linetype='dashed')+
            geom_line(aes(x=by,y=med),linewidth=1.3)+
            geom_point(aes(colour = by),size=4) +
            scale_colour_viridis_c(name='Year')+
            ggtitle(paste(title))+
            xlab("Year") + 
            ylab("log(Alpha)")+
            theme_classic(14)+
            theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                  strip.text = element_text(face="bold", size=12),
                  axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
          
          
        }
        
        plot1=ggplot2::ggplot(df, aes(S, R)) +
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = by_q[1]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = by_q[2]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,4],colour = by_q[3]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,5],colour = by_q[4]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,6],colour = by_q[5]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,7],colour = by_q[6]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = by_q[1]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = by_q[2]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,4],colour = by_q[3]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,5],colour = by_q[4]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,6],colour = by_q[5]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,7],colour = by_q[6]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,8],colour = by_q[7]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,9],colour = by_q[8]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,10],colour = by_q[9]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,11],colour = by_q[10]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,12],colour = by_q[11]),linewidth=1.3)+
          geom_point(aes(colour = by),size=2.5) +
          scale_colour_viridis_c(name='Year')+
          ggtitle(title)+
          xlab("Spawners") + 
          ylab("Recruits")+
          xlim(0, max(df$S))+
          ylim(0, max(df$R))+
          theme_classic(14)+
          theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                strip.text = element_text(face="bold", size=12),
                axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))

       
        legend = cowplot::get_legend(plot1)
        
        plot_rw_a=cowplot::plot_grid(plot1 + theme(legend.position="none"),
                      plot2 + theme(legend.position="none"),
                       ncol=2,nrow=1,labels=c("A","B"))
        plot=cowplot::plot_grid(plot_rw_a,legend,rel_widths = c(3,.25))
        if(sr_only==TRUE){plot=plot1}
      }
      if(par=='b'){ ###rw beta=====
          if(form=='stan'){
            post=rstan::extract(mod)
            pred_df=data.frame(x_new)
            
            for(n in 1:length(by_q)){
              pred_df[,1+n]=exp(median(post$log_a)-median(post$b[,match(by_q[n],df$by)])*x_new)*x_new
            }
            
            beta_df=data.frame(by=seq(min(df$by),max(df$by)),med=apply(post$S_max,2,median),l90=apply(post$S_max,2,quantile,0.1),u90=apply(post$S_max,2,quantile,0.9))
            
            plot2=ggplot2::ggplot(beta_df, aes(by,med)) +
              geom_line(aes(x=by,y=med),linewidth=1.3)+
              geom_point(aes(colour = by),size=4) +
              scale_colour_viridis_c(name='Year')+
              geom_ribbon(aes(ymin =l90, ymax =u90), alpha = 0.2)+
              ggtitle(paste(title))+
              xlab("Year") + 
              ylab("Smax")+
              theme_classic(14)+
              theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                    strip.text = element_text(face="bold", size=12),
                    axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
            
          }
          if(form=='tmb'){
            pred_df=data.frame(x_new)
            for(n in 1:length(by_q)){
              pred_df[,1+n]=exp(mod$alpha-mod$beta[match(by_q[n],df$by)]*x_new)*x_new
            }
            beta_df=data.frame(by=df$by,med=mod$Smax)
            
            plot2=ggplot2::ggplot(beta_df, aes(by,med)) +
             geom_line(aes(x=by,y=med),linewidth=1.3)+
              geom_point(aes(colour = by),size=4) +
              scale_colour_viridis_c(name='Year')+
              ggtitle(paste(title))+
              xlab("Year") + 
              ylab("Smax")+
              theme_classic(14)+
              theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                    strip.text = element_text(face="bold", size=12),
                    axis.text=element_text(face="bold",size=14),axis.title = element_text(face="bold",size=14),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
            
            
          }
          
        
        plot1=ggplot2::ggplot(df, aes(S, R)) +
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = by_q[1]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = by_q[2]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,4],colour = by_q[3]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,5],colour = by_q[4]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,6],colour = by_q[5]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,7],colour = by_q[6]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = by_q[1]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = by_q[2]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,4],colour = by_q[3]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,5],colour = by_q[4]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,6],colour = by_q[5]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,7],colour = by_q[6]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,8],colour = by_q[7]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,9],colour = by_q[8]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,10],colour = by_q[9]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,11],colour = by_q[10]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,12],colour = by_q[11]),linewidth=1.3)+
          geom_point(aes(colour = by),size=2.5) +
          scale_colour_viridis_c(name='Year')+
          ggtitle(title)+
          xlab("Spawners") + 
          ylab("Recruits")+
          xlim(0, max(df$S))+
          ylim(0, max(df$R))+
          theme_classic(14)+
          theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                strip.text = element_text(face="bold", size=12),
                axis.text=element_text(face="bold",size=14),axis.title = element_text(face="bold",size=14),plot.title = element_text(face = "bold", hjust = 0.5,size=15))

        
        legend = cowplot::get_legend(plot1)
        
        plot_rw_b=cowplot::plot_grid(plot1 + theme(legend.position="none"),
                                     plot2 + theme(legend.position="none"),
                                     ncol=2,nrow=1,labels=c("A","B"))
        plot=cowplot::plot_grid(plot_rw_b,legend,rel_widths = c(3,.25))
        if(sr_only==TRUE){plot=plot1}
      }
      if(par=='both'){ #rw alpha beta=====
        if(form=='stan'){
          post=rstan::extract(mod)
          pred_df=data.frame(x_new)
          for(n in 1:length(by_q)){
            pred_df[,1+n]=exp(median(post$log_a[,match(by_q[n],df$by)])-median(post$b[,match(by_q[n],df$by)])*x_new)*x_new
          }
          
          alphabeta_df=data.frame(by=seq(min(df$by),max(df$by)),a_med=apply(post$log_a,2,median),a_l90=apply(post$log_a,2,quantile,0.15),a_u90=apply(post$log_a,2,quantile,0.85),b_med=apply(post$S_max,2,median),b_l90=apply(post$S_max,2,quantile,0.1),b_u90=apply(post$S_max,2,quantile,0.9))
         
          plot2=ggplot2::ggplot(alphabeta_df, aes(by,a_med)) +
            geom_line(aes(x=by,y=a_med),linewidth=1.3)+
            geom_point(aes(colour = by),size=3) +
            scale_colour_viridis_c(name='Year')+
            geom_ribbon(aes(ymin =a_l90, ymax =a_u90), alpha = 0.2)+
            ggtitle(paste(title))+
            xlab("Year") + 
            ylab("log(Alpha)")+
            theme_classic(14)+
            theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                  strip.text = element_text(face="bold", size=12),
                  axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
          
          plot3=ggplot2::ggplot(alphabeta_df, aes(by,b_med)) +
            geom_line(aes(x=by,y=b_med),linewidth=1.3)+
            geom_point(aes(colour = by),size=3) +
            scale_colour_viridis_c(name='Year')+
            geom_ribbon(aes(ymin =b_l90, ymax =b_u90), alpha = 0.2)+
            ggtitle(paste(title))+
            xlab("Year") + 
            ylab("Smax")+
            theme_classic(14)+
            theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                  strip.text = element_text(face="bold", size=12),
                  axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
          
        }
        if(form=='tmb'){
          pred_df=data.frame(x_new)
          for(n in 1:length(by_q)){
            pred_df[,1+n]=exp(mod$alpha[match(by_q[n],df$by)]-mod$beta[match(by_q[n],df$by)]*x_new)*x_new
          }
          alphabeta_df=data.frame(by=df$by,a_med=mod$alpha,b_med=mod$Smax)
          
          plot2=ggplot2::ggplot(alphabeta_df, aes(by,a_med)) +
            geom_line(aes(x=by,y=a_med),linewidth=1.3)+
            geom_point(aes(colour = by),size=3) +
            scale_colour_viridis_c(name='Year')+
            ggtitle(paste(title))+
            xlab("Year") + 
            ylab("log(Alpha)")+
            theme_classic(14)+
            theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                  strip.text = element_text(face="bold", size=12),
                  axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
          
          plot3=ggplot2::ggplot(alphabeta_df, aes(by,b_med)) +
            geom_line(aes(x=by,y=b_med),linewidth=1.3)+
            geom_point(aes(colour = by),size=3) +
            scale_colour_viridis_c(name='Year')+
            ggtitle(paste(title))+
            xlab("Year") + 
            ylab("Smax")+
            theme_classic(14)+
            theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                  strip.text = element_text(face="bold", size=12),
                  axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
          
          
        }
        
        plot1=ggplot2::ggplot(df, aes(S, R)) +
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = by_q[1]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = by_q[2]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,4],colour = by_q[3]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,5],colour = by_q[4]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,6],colour = by_q[5]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,7],colour = by_q[6]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = by_q[1]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = by_q[2]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,4],colour = by_q[3]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,5],colour = by_q[4]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,6],colour = by_q[5]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,7],colour = by_q[6]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,8],colour = by_q[7]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,9],colour = by_q[8]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,10],colour = by_q[9]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,11],colour = by_q[10]),linewidth=1.3)+
          geom_line(data=pred_df,aes(x=x_new,y=pred_df[,12],colour = by_q[11]),linewidth=1.3)+
          geom_point(aes(colour = by),size=2.5) +
          scale_colour_viridis_c(name='Year')+
          ggtitle(title)+
          xlab("Spawners") + 
          ylab("Recruits")+
          xlim(0, max(df$S))+
          ylim(0, max(df$R))+
          theme_classic(14)+
          theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                strip.text = element_text(face="bold", size=12),
                axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
        
        
        legend = cowplot::get_legend(plot1)
        
        plot=cowplot::plot_grid(plot1 + theme(legend.position="none"),
                                     plot2 + theme(legend.position="none"),
                                     plot3 + theme(legend.position="none"),
                                     legend,
                                     ncol=2,nrow=2,labels=c("A","B","C"))
        if(sr_only==TRUE){plot=plot1}
        
      }
      
    }
    if(type=='hmm'){ 
        
        x_new=seq(min(df$S),max(df$S),length.out=200)
        pred_df=data.frame(x_new)
        
        if(par=='a'){ #hmm alpha====
          
          if(form=='stan'){
          post=rstan::extract(mod)
          pred_df[,2]=exp(median(post$log_a[,1])-median(post$b)*x_new)*x_new
          pred_df[,3]=exp(median(post$log_a[,2])-median(post$b)*x_new)*x_new
          df$gamma=apply(post$gamma[,,2],2,median)
          gamma_df=data.frame(by=df$by,gamma=apply(post$gamma[,,2],2,median),gamma_l90=apply(post$gamma[,,2],2,quantile,0.1),gamma_u90=apply(post$gamma[,,2],2,quantile,0.9))
          
          plot2=ggplot2::ggplot(gamma_df, aes(by,gamma)) +
            ylim(0,1)+
            geom_hline(yintercept=0.5,linetype='dashed')+
            geom_line(aes(x=by,y=gamma),linewidth=1.3)+
            geom_point(aes(colour = gamma),size=4) +
            scale_colour_viridis_c(name='p')+
            geom_ribbon(aes(ymin =gamma_l90, ymax =gamma_u90), alpha = 0.2)+
            ggtitle(paste(title))+
            xlab("Year") + 
            ylab("Prob. of high prod. regime")+
            theme_classic(14)+
            theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                  strip.text = element_text(face="bold", size=12),
                  axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
          }
          
          if(form=='tmb'){
            pred_df[,2]=exp(mod$alpha[1]-mod$beta*x_new)*x_new
            pred_df[,3]=exp(mod$alpha[2]-mod$beta*x_new)*x_new
            gamma_df=data.frame(by=df$by,gamma=mod$probregime[2,])
  
            plot2=ggplot2::ggplot(gamma_df, aes(by,gamma)) +
              ylim(0,1)+
              geom_hline(yintercept=0.5,linetype='dashed')+
              geom_line(aes(x=by,y=gamma),linewidth=1.3)+
              geom_point(aes(colour = gamma),size=4) +
              scale_colour_viridis_c(name='P(High prod. regime)')+
              ggtitle(paste(title))+
              xlab("Year") + 
              ylab("Prob. of high prod. regime")+
              theme_classic(14)+
              theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                    strip.text = element_text(face="bold", size=12),
                    axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
            
          }
          
          plot1=ggplot2::ggplot(df, aes(S, R)) +
            geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = min(gamma_df$gamma)),linewidth=1.3)+
            geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = max(gamma_df$gamma)),linewidth=1.3)+
            geom_point(aes(colour = gamma_df$gamma),size=2.5) +
            scale_colour_viridis_c(name='p')+
            ggtitle(title)+
            xlab("Spawners") + 
            ylab("Recruits")+
            theme_classic(14)+
            xlim(0, max(df$S))+
            ylim(0, max(df$R))+
            theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                  strip.text = element_text(face="bold", size=12),
                  axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
          
          
          
          legend = cowplot::get_legend(plot1)
          
          plot_hmm_a=cowplot::plot_grid(plot1 + theme(legend.position="none"),
                                       plot2 + theme(legend.position="none"),
                                       ncol=2,nrow=1,labels=c("A","B"))
          plot=cowplot::plot_grid(plot_hmm_a,legend,rel_widths = c(3,.3))
          if(sr_only==TRUE){plot=plot1}
        }
        
        if(par=='b'){ #hmm beta====
          
          if(form=='stan'){
          post=rstan::extract(mod)
          pred_df[,2]=exp(median(post$log_a)-median(post$b[,1])*x_new)*x_new
          pred_df[,3]=exp(median(post$log_a)-median(post$b[,2])*x_new)*x_new
          df$gamma=apply(post$gamma[,,2],2,median)
          gamma_df=data.frame(by=df$by,gamma=apply(post$gamma[,,1],2,median),gamma_l90=apply(post$gamma[,,1],2,quantile,0.1),gamma_u90=apply(post$gamma[,,1],2,quantile,0.9))
          
          plot2=ggplot2::ggplot(gamma_df, aes(by,gamma)) +
            ylim(0,1)+
            geom_hline(yintercept=0.5,linetype='dashed')+
            geom_line(aes(x=by,y=gamma),linewidth=1.3)+
            geom_point(aes(colour = gamma),size=4) +
            scale_colour_viridis_c(name='p')+
            geom_ribbon(aes(ymin =gamma_l90, ymax =gamma_u90), alpha = 0.2)+
            ggtitle(paste(title))+
            xlab("Year") + 
            ylab("Prob. of high cap. regime")+
            theme_classic(14)+
            theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                  strip.text = element_text(face="bold", size=12),
                  axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
          }
          if(form=='tmb'){
            pred_df[,2]=exp(mod$alpha-mod$beta[1]*x_new)*x_new
            pred_df[,3]=exp(mod$alpha-mod$beta[2]*x_new)*x_new
            gamma_df=data.frame(by=df$by,gamma=mod$probregime[1,])
            
            plot2=ggplot2::ggplot(gamma_df, aes(by,gamma)) +
              ylim(0,1)+
              geom_hline(yintercept=0.5,linetype='dashed')+
              geom_line(aes(x=by,y=gamma),linewidth=1.3)+
              geom_point(aes(colour = gamma),size=4) +
              scale_colour_viridis_c(name='p')+
              ggtitle(paste(title))+
              xlab("Year") + 
              ylab("Prob. of high cap. regime")+
              theme_classic(14)+
              theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                    strip.text = element_text(face="bold", size=12),
                    axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
            
          }
          
          plot1=ggplot2::ggplot(df, aes(S, R)) +
            geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = max(gamma_df$gamma)),linewidth=1.3)+
            geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = min(gamma_df$gamma)),linewidth=1.3)+
            geom_point(aes(colour = gamma_df$gamma),size=2.5) +
            scale_colour_viridis_c(name='p')+
            ggtitle(title)+
            xlab("Spawners") + 
            ylab("Recruits")+
            xlim(0, max(df$S))+
            ylim(0, max(df$R))+
            theme_classic(14)+
            theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                  strip.text = element_text(face="bold", size=12),
                  axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
          
         
          
          legend = cowplot::get_legend(plot1)
          
          plot_hmm_b=cowplot::plot_grid(plot1 + theme(legend.position="none"),
                                        plot2 + theme(legend.position="none"),
                                        ncol=2,nrow=1,labels=c("A","B"))
          plot=cowplot::plot_grid(plot_hmm_b,legend,rel_widths = c(3,.3))
          if(sr_only==TRUE){plot=plot1}
        }
        if(par=='both'){ #hmm alpha beta====
          if(form=='stan'){
          post=rstan::extract(mod)
          pred_df[,2]=exp(median(post$log_a[,1])-median(post$b[,1])*x_new)*x_new
          pred_df[,3]=exp(median(post$log_a[,2])-median(post$b[,2])*x_new)*x_new
          df$gamma=apply(post$gamma[,,2],2,median)
          gamma_df=data.frame(by=df$by,gamma=apply(post$gamma[,,2],2,median),gamma_l90=apply(post$gamma[,,2],2,quantile,0.1),gamma_u90=apply(post$gamma[,,2],2,quantile,0.9))
          
          plot2=ggplot2::ggplot(gamma_df, aes(by,gamma)) +
            ylim(0,1)+
            geom_hline(yintercept=0.5,linetype='dashed')+
            geom_line(aes(x=by,y=gamma),linewidth=1.3)+
            geom_point(aes(colour = gamma),size=4) +
            scale_colour_viridis_c(name='p')+
            geom_ribbon(aes(ymin =gamma_l90, ymax =gamma_u90), alpha = 0.2)+
            ggtitle(paste(title))+
            xlab("Year") + 
            ylab("Prob. of high prod. regime")+
            theme_classic(14)+
            theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                  strip.text = element_text(face="bold", size=12),
                  axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
          }
          if(form=='tmb'){
            pred_df[,2]=exp(mod$alpha[1]-mod$beta[1]*x_new)*x_new
            pred_df[,3]=exp(mod$alpha[2]-mod$beta[2]*x_new)*x_new
            gamma_df=data.frame(by=df$by,gamma=mod$probregime[1,])
            
            plot2=ggplot2::ggplot(gamma_df, aes(by,gamma)) +
              ylim(0,1)+
              geom_hline(yintercept=0.5,linetype='dashed')+
              geom_line(aes(x=by,y=gamma),linewidth=1.3)+
              geom_point(aes(colour = gamma),size=4) +
              scale_colour_viridis_c(name='p')+
              ggtitle(paste(title))+
              xlab("Year") + 
              ylab("Prob. of high cap. regime")+
              theme_classic(14)+
              theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                    strip.text = element_text(face="bold", size=12),
                    axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
          }
          
          plot1=ggplot2::ggplot(df, aes(S, R)) +
            geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = min(gamma_df$gamma)),linewidth=1.3)+
            geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = max(gamma_df$gamma)),linewidth=1.3)+
            geom_point(aes(colour = gamma_df$gamma),size=2.5) +
            scale_colour_viridis_c(name='p')+
            ggtitle(title)+
            xlab("Spawners") + 
            ylab("Recruits")+
            xlim(0, max(df$S))+
            ylim(0, max(df$R))+
            theme_classic(14)+
            theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                  strip.text = element_text(face="bold", size=12),
                  axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
          
          legend = cowplot::get_legend(plot1)
          
          plot_hmm_ab=cowplot::plot_grid(plot1 + theme(legend.position="none"),
                                        plot2 + theme(legend.position="none"),
                                        ncol=2,nrow=1,labels=c("A","B"))
          plot=cowplot::plot_grid(plot_hmm_ab,legend,rel_widths = c(3,.3))
          if(sr_only==TRUE){plot=plot1}
        }
    }
    if(make.pdf==TRUE){
      if(type=='static'&ac==FALSE){ pdf(here(path,paste(paste(title,type,form,sep='_'),'.pdf',sep='')),width=8,height=6)}
      if(type=='static'&ac==TRUE){ pdf(here(path,paste(paste(title,type,'ac',form,sep='_'),'.pdf',sep='')),width=8,height=6)}
      if(type=='rw'&par=='both'){ pdf(here(path,paste(paste(title,type,par,form,sep='_'),'.pdf',sep='')),width=10,height=10)}
      if(type=='rw'&par!='both'){ pdf(here(path,paste(paste(title,type,par,form,sep='_'),'.pdf',sep='')),width=14,height=6)}
      if(type=='hmm'){ pdf(here(path,paste(paste(title,type,par,form,sep='_'),'.pdf',sep='')),width=14,height=6)}
      return(plot)
      
      dev.off()
      dev.off()
    }
  if(make.pdf==FALSE){
    return(plot)
  }
  
  }

#' Prior predictive check 
#'
#' Any arguments to pass to [rstan::sampling].
#'
#' @param fit A model fit in rstan from e.g. ricker_stan, ricker_rw_stan, ricker_hmm_stan
#'
#' @export
prior_check<- function(data,type=c('static','rw'),AC=FALSE,par=c('a','b','both'),control = stancontrol()){
  if(type=='static'&AC==FALSE){
   
    sm=sr_prior_mod(type='static')
      
   if(is.null(smax_priors)==T){
      datm = list(N=nrow(data),
                  L=max(data$by)-min(data$by)+1,
                  ii=seq_len(nrow(data)),
                  R_S =data$logRS,
                  S=data$S,
                  pSmax_mean=max(data$S)/2,
                  pSmax_sig=max(data$S))
      
    }
    if(is.null(smax_priors)==F){
      datm = list(N=nrow(data),
                  L=max(data$by)-min(data$by)+1,
                  ii=seq_len(nrow(data)),
                  R_S =data$logRS,
                  S=data$S,
                  pSmax_mean=smax_priors[1],
                  pSmax_sig=smax_priors[2])
    }
  }
  if(type=='static'&AC==TRUE){
    if(is.null(smax_priors)==T){
      sm=sr_prior_mod(type='static',AC=TRUE)
    }
    if(is.null(smax_priors)==F){
      sm=sr_prior_mod(type='static',AC=TRUE,smax_priors=smax_priors)
    }
  }
  
  if(type=='rw'&par=='a'){
    if(is.null(smax_priors)==T){
      sm=sr_prior_mod(type='rw',par='a')
    }
    if(is.null(smax_priors)==F){
      sm=sr_prior_mod(type='rw',par='a',smax_priors=smax_priors)
    }
  }
  if(type=='rw'&par=='b'){
    if(is.null(smax_priors)==T){
      sm=sr_prior_mod(type='rw',par='b')
    }
    if(is.null(smax_priors)==F){
      sm=sr_prior_mod(type='rw',par='b',smax_priors=smax_priors)
    }
  }
  if(type=='rw'&par=='both'){
    if(is.null(smax_priors)==T){
      sm=sr_prior_mod(type='rw',par='both')
    }
    if(is.null(smax_priors)==F){
      sm=sr_prior_mod(type='rw',par='both',smax_priors=smax_priors)
    }
  }
  
  fit<-rstan::sampling(sm, data=datm,
                       control = control, warmup = warmup, 
                       chains = chains, iter = iter,verbose=FALSE)
  
  yrep=rstan::extract(fit,pars='y_rep')
  
  
  
}


#' Posterior predictive check plot
#'
#' @param fit A model fit in rstan from e.g. ricker_stan, ricker_rw_stan, ricker_hmm_stan
#' @param data the data that were used for the model in 'fit'
#' @export
post_check<- function(fit,data,pdf=FALSE,path=NULL,filename=NULL){
  yrep_RS=rstan::extract(fit$summary,pars='y_rep',permuted=FALSE)
  yrep_R=array(NA,dim=dim(yrep_RS))
  for(t in 1:6){
    yrep_R[,t,]=log10(exp(yrep_RS[,t,])*data$S)
  }
  
  emp_RS=density(data$logRS,bw=0.02) #empirical distribution of log(R/S)
  emp_R=density(log10(data$R),bw=0.02) #empirical distribution of recruits
  
  cols=RColorBrewer::brewer.pal(n=7,'Blues')
  
  if(pdf==T){pdf(paste(path,'/',filename,'.pdf',sep=''),height=12,width=8)}
  par(mfrow=c(2,2))
  
  hist(data$logRS,main='',xaxt='n',breaks=30,freq=T,xlab='',ylab='counts of observations in time-series',col=adjustcolor('darkred',alpha.f=0.6),border='white',xlim=c(min(yrep_RS),max(yrep_RS)),xaxt='n')
  text('empirical obs.',x=par('usr')[2]-((par('usr')[2]-par('usr')[1])*0.2),y=par('usr')[4]-((par('usr')[4]-par('usr')[3])*0.1),col=adjustcolor('darkred',alpha.f=0.6))
  par(new=T)
  plot(rep(0,length(emp_RS$x))~emp_RS$x,type='n',col='darkred',ylab='',xlab='log(R/S)',xlim=c(min(yrep_RS),max(yrep_RS)),ylim=c(0,1),yaxt='n',bty='l')
   for(i in 1:6){
    d=density(yrep_RS[,i,],bw=0.02)
    lines(d$y/max(d$y)~d$x,col=cols[i+1])
    text(paste('pred. chain',i,sep=' '),x=par('usr')[2]-((par('usr')[2]-par('usr')[1])*0.2),y=par('usr')[4]-((par('usr')[4]-par('usr')[3])*(0.1+0.05*i)),col=cols[i+1])
  }
  hist(log10(data$R),main='',xaxt='n',breaks=30,freq=T,xlab='',ylab='counts of observations in time-series',col=adjustcolor('darkred',alpha.f=0.6),border='white',xlim=c(min(yrep_R),max(yrep_R)),xaxt='n')
  text('empirical obs.',x=par('usr')[2]-((par('usr')[2]-par('usr')[1])*0.2),y=par('usr')[4]-((par('usr')[4]-par('usr')[3])*0.1),col=adjustcolor('darkred',alpha.f=0.6))
  par(new=T)
  plot(rep(0,length(emp_R$x))~emp_R$x,type='n',col='darkred',ylab='',xlab='recuits (log10 axis)',xlim=c(min(yrep_R),max(yrep_R)),ylim=c(0,1),yaxt='n',xaxt='n',bty='l')
  for(i in 1:6){
    d=density(yrep_R[,i,],bw=0.02)
    lines(d$y/max(d$y)~d$x,col=cols[i+1])
    text(paste('pred. chain',i,sep=' '),x=par('usr')[2]-((par('usr')[2]-par('usr')[1])*0.2),y=par('usr')[4]-((par('usr')[4]-par('usr')[3])*(0.1+0.05*i)),col=cols[i+1])
  }
  pow <- c(min(round(log10(data$R)))-1):c(max(round(log10(data$R))))
  axis(1, col="black", at=seq(min(round(log10(data$R)))-1,max(round(log10(data$R))),by=1),tcl=-0.45, cex.axis=1.2,labels=10^pow)
  ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
  axis(1, log10(ticksat), col="black", labels=NA,
       tcl=-0.2, lwd=0, lwd.ticks=1)
  
  smaxs=rstan::extract(fit$summary,pars=c('prior_Smax','S_max'))
  hist(c(smaxs$prior_Smax/1e3),breaks=30,freq=T,xlim=c(0,max(c(smaxs[[1]]/1e3,smaxs[[2]]/1e3))),xlab='spawners (1000s)',col=adjustcolor('darkorange',alpha.f=0.5),border='white',main='')
  par(new=T)
  hist(c(smaxs$S_max/1e3),breaks=30,freq=T,xlim=c(0,max(c(smaxs[[1]]/1e3,smaxs[[2]]/1e3))),xlab='spawners (1000s)',col=adjustcolor('navy',alpha.f=0.5),border='white',main='',yaxt='n')
  text('Smax prior',x=par('usr')[2]-((par('usr')[2]-par('usr')[1])*0.2),y=par('usr')[4]-((par('usr')[4]-par('usr')[3])*0.1),col='darkorange')
  text('Smax posterior',x=par('usr')[2]-((par('usr')[2]-par('usr')[1])*0.2),y=par('usr')[4]-((par('usr')[4]-par('usr')[3])*0.2),col='navy')
  
  plot(c(data$R/1e3)~c(data$S/1e3),xlab='spawners (1000s)',ylab='recruits (1000s)',type='n',ylim=c(0,max(data$R/1e3)),xlim=c(0,max(data$S/1e3)*1.2),bty='l')
  sn=seq(0,max(data$S))
  muR=exp(fit$alpha[1]-fit$beta[1]*sn)*sn
  lines(c(muR/1e3)~c(sn/1e3))
  par(new=T)
  hist(c(smaxs$prior_Smax/1e3),breaks=30,freq=T,xlim=c(0,max(data$S/1e3)*1.2),ylim=c(0,1e3),xlab='',col=adjustcolor('darkorange',alpha.f=0.5),border='white',main='',xaxt='n',yaxt='n',ylab='')
  par(new=T)
  hist(c(smaxs$S_max/1e3),breaks=30,freq=T,xlim=c(0,max(data$S/1e3)*1.2),ylim=c(0,1e3),xlab='',col=adjustcolor('navy',alpha.f=0.5),border='white',main='',yaxt='n',xaxt='n',ylab='')
  text(y=c(data$R/1e3),x=c(data$S/1e3),data$by,col='darkred',cex=0.8,font=2) 
  if(pdf==T){
    dev.off()
    dev.off()
    }
}

sr_plot2=function(df,mod,title,make.pdf=FALSE,path,type=c('static','rw','hmm'),par=c('a','b','both'),form=c('stan','tmb'),ac=FALSE,sr_only=FALSE){
  if(type=='static'){ #static====
    x_new=seq(min(df$S),max(df$S),length.out=200)
    
    if(form=='stan'){
      post=as.data.frame(mod$draws(format='draws_matrix'))
      pred_df=data.frame(pred=exp(median(post$log_a)-median(post$b)*x_new)*x_new,x_new=x_new)
    }
    if(form=='tmb'){
      pred_df=data.frame(pred=exp(mod$alpha-mod$beta*x_new)*x_new,x_new=x_new)
    }
    plot=ggplot2::ggplot(df, aes(S, R)) +
      geom_line(data=pred_df,aes(x=x_new,y=pred),linewidth=1.3)+
      geom_point(aes(colour = by),size=2.5) +
      scale_colour_viridis_c(name='Year')+
      ggtitle(title)+
      xlab("Spawners") + 
      ylab("Recruits")+
      xlim(0, max(df$S))+
      ylim(0, max(df$R))+
      theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
            strip.text = element_text(face="bold", size=12),
            axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
  }
  if(type=='rw'){
    x_new=seq(min(df$S),max(df$S),length.out=200)
    by_q=round(quantile(df$by,seq(0,1,by=0.1)))
    if(par=='a'){ #rw alpha=====
      if(form=='stan'){
        log_a=as.data.frame(mod$draws(variables='log_a',format='draws_matrix'))
        b=as.data.frame(mod$draws(variables='b',format='draws_matrix'))
        pred_df=data.frame(x_new)
        for(n in 1:length(by_q)){
          pred_df[,1+n]=exp(median(log_a[,match(by_q[n],df$by)])-median(b$b)*x_new)*x_new
        }
        alpha_df=data.frame(by=seq(min(df$by),max(df$by)),med=apply(log_a,2,median),l90=apply(log_a,2,quantile,0.1),u90=apply(log_a,2,quantile,0.9))
        plot2=ggplot2::ggplot(alpha_df, aes(by,med)) +
          geom_line(aes(x=by,y=med),linewidth=1.3)+
          geom_point(aes(colour = by),size=4) +
          scale_colour_viridis_c(name='Year')+
          geom_ribbon(aes(ymin =l90, ymax =u90), alpha = 0.2)+
          ggtitle(title)+
          xlab("Year") + 
          ylab("log(Alpha)")+
          theme_classic(14)+
          theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                strip.text = element_text(face="bold", size=12),
                axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
        
      }
      
      if(form=='tmb'){
        pred_df=data.frame(x_new)
        for(n in 1:length(by_q)){
          pred_df[,1+n]=exp(mod$alpha[match(by_q[n],df$by)]-mod$beta*x_new)*x_new
        }
        alpha_df=data.frame(by=df$by,med=mod$alpha)
        
        plot2=ggplot2::ggplot(alpha_df, aes(by,med)) +
          geom_hline(yintercept=0,linetype='dashed')+
          geom_line(aes(x=by,y=med),linewidth=1.3)+
          geom_point(aes(colour = by),size=4) +
          scale_colour_viridis_c(name='Year')+
          ggtitle(paste(title))+
          xlab("Year") + 
          ylab("log(Alpha)")+
          theme_classic(14)+
          theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                strip.text = element_text(face="bold", size=12),
                axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
        
        
      }
      
      plot1=ggplot2::ggplot(df, aes(S, R)) +
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = by_q[1]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = by_q[2]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,4],colour = by_q[3]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,5],colour = by_q[4]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,6],colour = by_q[5]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,7],colour = by_q[6]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = by_q[1]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = by_q[2]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,4],colour = by_q[3]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,5],colour = by_q[4]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,6],colour = by_q[5]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,7],colour = by_q[6]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,8],colour = by_q[7]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,9],colour = by_q[8]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,10],colour = by_q[9]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,11],colour = by_q[10]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,12],colour = by_q[11]),linewidth=1.3)+
        geom_point(aes(colour = by),size=2.5) +
        scale_colour_viridis_c(name='Year')+
        ggtitle(title)+
        xlab("Spawners") + 
        ylab("Recruits")+
        xlim(0, max(df$S))+
        ylim(0, max(df$R))+
        theme_classic(14)+
        theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
              strip.text = element_text(face="bold", size=12),
              axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
      
      
      legend = cowplot::get_legend(plot1)
      
      plot_rw_a=cowplot::plot_grid(plot1 + theme(legend.position="none"),
                                   plot2 + theme(legend.position="none"),
                                   ncol=2,nrow=1,labels=c("A","B"))
      plot=cowplot::plot_grid(plot_rw_a,legend,rel_widths = c(3,.25))
      if(sr_only==TRUE){plot=plot1}
    }
    if(par=='b'){ ###rw beta=====
      if(form=='stan'){
        post=as.data.frame(mod$draws(format='draws_matrix'))
        pred_df=data.frame(x_new)
        
        for(n in 1:length(by_q)){
          pred_df[,1+n]=exp(median(post$log_a)-median(post$b[,match(by_q[n],df$by)])*x_new)*x_new
        }
        
        beta_df=data.frame(by=df$by,med=apply(post$S_max,2,median),l90=apply(post$S_max,2,quantile,0.1),u90=apply(post$S_max,2,quantile,0.9))
        
        plot2=ggplot2::ggplot(beta_df, aes(by,med)) +
          geom_line(aes(x=by,y=med),linewidth=1.3)+
          geom_point(aes(colour = by),size=4) +
          scale_colour_viridis_c(name='Year')+
          geom_ribbon(aes(ymin =l90, ymax =u90), alpha = 0.2)+
          ggtitle(paste(title))+
          xlab("Year") + 
          ylab("Smax")+
          theme_classic(14)+
          theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                strip.text = element_text(face="bold", size=12),
                axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
        
      }
      if(form=='tmb'){
        pred_df=data.frame(x_new)
        for(n in 1:length(by_q)){
          pred_df[,1+n]=exp(mod$alpha-mod$beta[match(by_q[n],df$by)]*x_new)*x_new
        }
        beta_df=data.frame(by=df$by,med=mod$Smax)
        
        plot2=ggplot2::ggplot(beta_df, aes(by,med)) +
          geom_line(aes(x=by,y=med),linewidth=1.3)+
          geom_point(aes(colour = by),size=4) +
          scale_colour_viridis_c(name='Year')+
          ggtitle(paste(title))+
          xlab("Year") + 
          ylab("Smax")+
          theme_classic(14)+
          theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                strip.text = element_text(face="bold", size=12),
                axis.text=element_text(face="bold",size=14),axis.title = element_text(face="bold",size=14),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
        
        
      }
      
      
      plot1=ggplot2::ggplot(df, aes(S, R)) +
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = by_q[1]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = by_q[2]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,4],colour = by_q[3]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,5],colour = by_q[4]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,6],colour = by_q[5]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,7],colour = by_q[6]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = by_q[1]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = by_q[2]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,4],colour = by_q[3]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,5],colour = by_q[4]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,6],colour = by_q[5]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,7],colour = by_q[6]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,8],colour = by_q[7]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,9],colour = by_q[8]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,10],colour = by_q[9]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,11],colour = by_q[10]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,12],colour = by_q[11]),linewidth=1.3)+
        geom_point(aes(colour = by),size=2.5) +
        scale_colour_viridis_c(name='Year')+
        ggtitle(title)+
        xlab("Spawners") + 
        ylab("Recruits")+
        xlim(0, max(df$S))+
        ylim(0, max(df$R))+
        theme_classic(14)+
        theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
              strip.text = element_text(face="bold", size=12),
              axis.text=element_text(face="bold",size=14),axis.title = element_text(face="bold",size=14),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
      
      
      legend = cowplot::get_legend(plot1)
      
      plot_rw_b=cowplot::plot_grid(plot1 + theme(legend.position="none"),
                                   plot2 + theme(legend.position="none"),
                                   ncol=2,nrow=1,labels=c("A","B"))
      plot=cowplot::plot_grid(plot_rw_b,legend,rel_widths = c(3,.25))
      if(sr_only==TRUE){plot=plot1}
    }
    if(par=='both'){ #rw alpha beta=====
      if(form=='stan'){
        post=as.data.frame(mod$draws(format='draws_matrix'))
        pred_df=data.frame(x_new)
        for(n in 1:length(by_q)){
          pred_df[,1+n]=exp(median(post$log_a[,match(by_q[n],df$by)])-median(post$b[,match(by_q[n],df$by)])*x_new)*x_new
        }
        
        alphabeta_df=data.frame(by=df$by,a_med=apply(post$log_a,2,median),a_l90=apply(post$log_a,2,quantile,0.15),a_u90=apply(post$log_a,2,quantile,0.85),b_med=apply(post$S_max,2,median),b_l90=apply(post$S_max,2,quantile,0.1),b_u90=apply(post$S_max,2,quantile,0.9))
        
        plot2=ggplot2::ggplot(alphabeta_df, aes(by,a_med)) +
          geom_line(aes(x=by,y=a_med),linewidth=1.3)+
          geom_point(aes(colour = by),size=3) +
          scale_colour_viridis_c(name='Year')+
          geom_ribbon(aes(ymin =a_l90, ymax =a_u90), alpha = 0.2)+
          ggtitle(paste(title))+
          xlab("Year") + 
          ylab("log(Alpha)")+
          theme_classic(14)+
          theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                strip.text = element_text(face="bold", size=12),
                axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
        
        plot3=ggplot2::ggplot(alphabeta_df, aes(by,b_med)) +
          geom_line(aes(x=by,y=b_med),linewidth=1.3)+
          geom_point(aes(colour = by),size=3) +
          scale_colour_viridis_c(name='Year')+
          geom_ribbon(aes(ymin =b_l90, ymax =b_u90), alpha = 0.2)+
          ggtitle(paste(title))+
          xlab("Year") + 
          ylab("Smax")+
          theme_classic(14)+
          theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                strip.text = element_text(face="bold", size=12),
                axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
        
      }
      if(form=='tmb'){
        pred_df=data.frame(x_new)
        for(n in 1:length(by_q)){
          pred_df[,1+n]=exp(mod$alpha[match(by_q[n],df$by)]-mod$beta[match(by_q[n],df$by)]*x_new)*x_new
        }
        alphabeta_df=data.frame(by=df$by,a_med=mod$alpha,b_med=mod$Smax)
        
        plot2=ggplot2::ggplot(alphabeta_df, aes(by,a_med)) +
          geom_line(aes(x=by,y=a_med),linewidth=1.3)+
          geom_point(aes(colour = by),size=3) +
          scale_colour_viridis_c(name='Year')+
          ggtitle(paste(title))+
          xlab("Year") + 
          ylab("log(Alpha)")+
          theme_classic(14)+
          theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                strip.text = element_text(face="bold", size=12),
                axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
        
        plot3=ggplot2::ggplot(alphabeta_df, aes(by,b_med)) +
          geom_line(aes(x=by,y=b_med),linewidth=1.3)+
          geom_point(aes(colour = by),size=3) +
          scale_colour_viridis_c(name='Year')+
          ggtitle(paste(title))+
          xlab("Year") + 
          ylab("Smax")+
          theme_classic(14)+
          theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                strip.text = element_text(face="bold", size=12),
                axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
        
        
      }
      
      plot1=ggplot2::ggplot(df, aes(S, R)) +
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = by_q[1]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = by_q[2]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,4],colour = by_q[3]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,5],colour = by_q[4]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,6],colour = by_q[5]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,7],colour = by_q[6]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = by_q[1]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = by_q[2]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,4],colour = by_q[3]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,5],colour = by_q[4]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,6],colour = by_q[5]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,7],colour = by_q[6]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,8],colour = by_q[7]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,9],colour = by_q[8]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,10],colour = by_q[9]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,11],colour = by_q[10]),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,12],colour = by_q[11]),linewidth=1.3)+
        geom_point(aes(colour = by),size=2.5) +
        scale_colour_viridis_c(name='Year')+
        ggtitle(title)+
        xlab("Spawners") + 
        ylab("Recruits")+
        xlim(0, max(df$S))+
        ylim(0, max(df$R))+
        theme_classic(14)+
        theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
              strip.text = element_text(face="bold", size=12),
              axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
      
      
      legend = cowplot::get_legend(plot1)
      
      plot=cowplot::plot_grid(plot1 + theme(legend.position="none"),
                              plot2 + theme(legend.position="none"),
                              plot3 + theme(legend.position="none"),
                              legend,
                              ncol=2,nrow=2,labels=c("A","B","C"))
      if(sr_only==TRUE){plot=plot1}
      
    }
    
  }
  if(type=='hmm'){ 
    
    x_new=seq(min(df$S),max(df$S),length.out=200)
    pred_df=data.frame(x_new)
    
    if(par=='a'){ #hmm alpha====
      
      if(form=='stan'){
        post=as.data.frame(mod$draws(format='draws_matrix'))
        pred_df[,2]=exp(median(post$log_a[,1])-median(post$b)*x_new)*x_new
        pred_df[,3]=exp(median(post$log_a[,2])-median(post$b)*x_new)*x_new
        df$gamma=apply(post$gamma[,,2],2,median)
        gamma_df=data.frame(by=df$by,gamma=apply(post$gamma[,,2],2,median),gamma_l90=apply(post$gamma[,,2],2,quantile,0.1),gamma_u90=apply(post$gamma[,,2],2,quantile,0.9))
        
        plot2=ggplot2::ggplot(gamma_df, aes(by,gamma)) +
          ylim(0,1)+
          geom_hline(yintercept=0.5,linetype='dashed')+
          geom_line(aes(x=by,y=gamma),linewidth=1.3)+
          geom_point(aes(colour = gamma),size=4) +
          scale_colour_viridis_c(name='p')+
          geom_ribbon(aes(ymin =gamma_l90, ymax =gamma_u90), alpha = 0.2)+
          ggtitle(paste(title))+
          xlab("Year") + 
          ylab("Prob. of high prod. regime")+
          theme_classic(14)+
          theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                strip.text = element_text(face="bold", size=12),
                axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
      }
      
      if(form=='tmb'){
        pred_df[,2]=exp(mod$alpha[1]-mod$beta*x_new)*x_new
        pred_df[,3]=exp(mod$alpha[2]-mod$beta*x_new)*x_new
        gamma_df=data.frame(by=df$by,gamma=mod$probregime[2,])
        
        plot2=ggplot2::ggplot(gamma_df, aes(by,gamma)) +
          ylim(0,1)+
          geom_hline(yintercept=0.5,linetype='dashed')+
          geom_line(aes(x=by,y=gamma),linewidth=1.3)+
          geom_point(aes(colour = gamma),size=4) +
          scale_colour_viridis_c(name='P(High prod. regime)')+
          ggtitle(paste(title))+
          xlab("Year") + 
          ylab("Prob. of high prod. regime")+
          theme_classic(14)+
          theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                strip.text = element_text(face="bold", size=12),
                axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
        
      }
      
      plot1=ggplot2::ggplot(df, aes(S, R)) +
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = min(gamma_df$gamma)),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = max(gamma_df$gamma)),linewidth=1.3)+
        geom_point(aes(colour = gamma_df$gamma),size=2.5) +
        scale_colour_viridis_c(name='p')+
        ggtitle(title)+
        xlab("Spawners") + 
        ylab("Recruits")+
        theme_classic(14)+
        xlim(0, max(df$S))+
        ylim(0, max(df$R))+
        theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
              strip.text = element_text(face="bold", size=12),
              axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
      
      
      
      legend = cowplot::get_legend(plot1)
      
      plot_hmm_a=cowplot::plot_grid(plot1 + theme(legend.position="none"),
                                    plot2 + theme(legend.position="none"),
                                    ncol=2,nrow=1,labels=c("A","B"))
      plot=cowplot::plot_grid(plot_hmm_a,legend,rel_widths = c(3,.3))
      if(sr_only==TRUE){plot=plot1}
    }
    
    if(par=='b'){ #hmm beta====
      
      if(form=='stan'){
        post=as.data.frame(mod$draws(format='draws_matrix'))
        pred_df[,2]=exp(median(post$log_a)-median(post$b[,1])*x_new)*x_new
        pred_df[,3]=exp(median(post$log_a)-median(post$b[,2])*x_new)*x_new
        df$gamma=apply(post$gamma[,,2],2,median)
        gamma_df=data.frame(by=df$by,gamma=apply(post$gamma[,,1],2,median),gamma_l90=apply(post$gamma[,,1],2,quantile,0.1),gamma_u90=apply(post$gamma[,,1],2,quantile,0.9))
        
        plot2=ggplot2::ggplot(gamma_df, aes(by,gamma)) +
          ylim(0,1)+
          geom_hline(yintercept=0.5,linetype='dashed')+
          geom_line(aes(x=by,y=gamma),linewidth=1.3)+
          geom_point(aes(colour = gamma),size=4) +
          scale_colour_viridis_c(name='p')+
          geom_ribbon(aes(ymin =gamma_l90, ymax =gamma_u90), alpha = 0.2)+
          ggtitle(paste(title))+
          xlab("Year") + 
          ylab("Prob. of high cap. regime")+
          theme_classic(14)+
          theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                strip.text = element_text(face="bold", size=12),
                axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
      }
      if(form=='tmb'){
        pred_df[,2]=exp(mod$alpha-mod$beta[1]*x_new)*x_new
        pred_df[,3]=exp(mod$alpha-mod$beta[2]*x_new)*x_new
        gamma_df=data.frame(by=df$by,gamma=mod$probregime[1,])
        
        plot2=ggplot2::ggplot(gamma_df, aes(by,gamma)) +
          ylim(0,1)+
          geom_hline(yintercept=0.5,linetype='dashed')+
          geom_line(aes(x=by,y=gamma),linewidth=1.3)+
          geom_point(aes(colour = gamma),size=4) +
          scale_colour_viridis_c(name='p')+
          ggtitle(paste(title))+
          xlab("Year") + 
          ylab("Prob. of high cap. regime")+
          theme_classic(14)+
          theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                strip.text = element_text(face="bold", size=12),
                axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
        
      }
      
      plot1=ggplot2::ggplot(df, aes(S, R)) +
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = max(gamma_df$gamma)),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = min(gamma_df$gamma)),linewidth=1.3)+
        geom_point(aes(colour = gamma_df$gamma),size=2.5) +
        scale_colour_viridis_c(name='p')+
        ggtitle(title)+
        xlab("Spawners") + 
        ylab("Recruits")+
        xlim(0, max(df$S))+
        ylim(0, max(df$R))+
        theme_classic(14)+
        theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
              strip.text = element_text(face="bold", size=12),
              axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
      
      
      
      legend = cowplot::get_legend(plot1)
      
      plot_hmm_b=cowplot::plot_grid(plot1 + theme(legend.position="none"),
                                    plot2 + theme(legend.position="none"),
                                    ncol=2,nrow=1,labels=c("A","B"))
      plot=cowplot::plot_grid(plot_hmm_b,legend,rel_widths = c(3,.3))
      if(sr_only==TRUE){plot=plot1}
    }
    if(par=='both'){ #hmm alpha beta====
      if(form=='stan'){
        post=as.data.frame(mod$draws(format='draws_matrix'))
        pred_df[,2]=exp(median(post$log_a[,1])-median(post$b[,1])*x_new)*x_new
        pred_df[,3]=exp(median(post$log_a[,2])-median(post$b[,2])*x_new)*x_new
        df$gamma=apply(post$gamma[,,2],2,median)
        gamma_df=data.frame(by=df$by,gamma=apply(post$gamma[,,2],2,median),gamma_l90=apply(post$gamma[,,2],2,quantile,0.1),gamma_u90=apply(post$gamma[,,2],2,quantile,0.9))
        
        plot2=ggplot2::ggplot(gamma_df, aes(by,gamma)) +
          ylim(0,1)+
          geom_hline(yintercept=0.5,linetype='dashed')+
          geom_line(aes(x=by,y=gamma),linewidth=1.3)+
          geom_point(aes(colour = gamma),size=4) +
          scale_colour_viridis_c(name='p')+
          geom_ribbon(aes(ymin =gamma_l90, ymax =gamma_u90), alpha = 0.2)+
          ggtitle(paste(title))+
          xlab("Year") + 
          ylab("Prob. of high prod. regime")+
          theme_classic(14)+
          theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                strip.text = element_text(face="bold", size=12),
                axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
      }
      if(form=='tmb'){
        pred_df[,2]=exp(mod$alpha[1]-mod$beta[1]*x_new)*x_new
        pred_df[,3]=exp(mod$alpha[2]-mod$beta[2]*x_new)*x_new
        gamma_df=data.frame(by=df$by,gamma=mod$probregime[1,])
        
        plot2=ggplot2::ggplot(gamma_df, aes(by,gamma)) +
          ylim(0,1)+
          geom_hline(yintercept=0.5,linetype='dashed')+
          geom_line(aes(x=by,y=gamma),linewidth=1.3)+
          geom_point(aes(colour = gamma),size=4) +
          scale_colour_viridis_c(name='p')+
          ggtitle(paste(title))+
          xlab("Year") + 
          ylab("Prob. of high cap. regime")+
          theme_classic(14)+
          theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                strip.text = element_text(face="bold", size=12),
                axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
      }
      
      plot1=ggplot2::ggplot(df, aes(S, R)) +
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2],colour = min(gamma_df$gamma)),linewidth=1.3)+
        geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3],colour = max(gamma_df$gamma)),linewidth=1.3)+
        geom_point(aes(colour = gamma_df$gamma),size=2.5) +
        scale_colour_viridis_c(name='p')+
        ggtitle(title)+
        xlab("Spawners") + 
        ylab("Recruits")+
        xlim(0, max(df$S))+
        ylim(0, max(df$R))+
        theme_classic(14)+
        theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
              strip.text = element_text(face="bold", size=12),
              axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
      
      legend = cowplot::get_legend(plot1)
      
      plot_hmm_ab=cowplot::plot_grid(plot1 + theme(legend.position="none"),
                                     plot2 + theme(legend.position="none"),
                                     ncol=2,nrow=1,labels=c("A","B"))
      plot=cowplot::plot_grid(plot_hmm_ab,legend,rel_widths = c(3,.3))
      if(sr_only==TRUE){plot=plot1}
    }
  }
  if(make.pdf==TRUE){
    if(type=='static'&ac==FALSE){ pdf(here(path,paste(paste(title,type,form,sep='_'),'.pdf',sep='')),width=8,height=6)}
    if(type=='static'&ac==TRUE){ pdf(here(path,paste(paste(title,type,'ac',form,sep='_'),'.pdf',sep='')),width=8,height=6)}
    if(type=='rw'&par=='both'){ pdf(here(path,paste(paste(title,type,par,form,sep='_'),'.pdf',sep='')),width=10,height=10)}
    if(type=='rw'&par!='both'){ pdf(here(path,paste(paste(title,type,par,form,sep='_'),'.pdf',sep='')),width=14,height=6)}
    if(type=='hmm'&par=='both'){ pdf(here(path,paste(paste(title,type,par,form,sep='_'),'.pdf',sep='')),width=8,height=6)}
    if(type=='hmm'&par!='both'){ pdf(here(path,paste(paste(title,type,par,form,sep='_'),'.pdf',sep='')),width=14,height=6)}
    print(plot)
    dev.off()
  }
  if(make.pdf==FALSE){
    return(plot)
  }
  
}
