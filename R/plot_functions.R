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

sr_plot=function(df,mod,title,make.pdf=FALSE,path,type=c('static','rw','hmm'),par=c('a','b','both'),form=c('stan','tmb')){
  if(type=='static'){ #static====
    x_new=seq(min(df$S),max(df$S),length.out=200)

    if(form=='stan'){
      post=rstan::extract(mod)
      pred_df=data.frame(pred=exp(median(post$log_a)-median(post$b)*x_new)*x_new)
      }
    if(form=='tmb'){
      pred_df=data.frame(pred=exp(mod$alpha-mod$beta*x_new)*x_new)
    }
    plot=ggplot2::ggplot(df, aes(S, R)) +
      geom_line(data=pred_df,aes(x=x_new,y=pred),linewidth=1.3)+
      geom_point(aes(colour = by),size=2.5) +
      scale_colour_viridis_c(name='Year')+
      xlab("Spawners") + 
      ylab("Recruits")+
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
          alpha_df=data.frame(by=df$by,med=apply(post$log_a,2,median),l90=apply(post$log_a,2,quantile,0.1),u90=apply(post$log_a,2,quantile,0.9))
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
          theme_classic(14)+
          theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                strip.text = element_text(face="bold", size=12),
                axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))

       
        legend = cowplot::get_legend(plot1)
        
        plot_rw_a=cowplot::plot_grid(plot1 + theme(legend.position="none"),
                      plot2 + theme(legend.position="none"),
                       ncol=2,nrow=1,labels=c("A","B"))
        plot=cowplot::plot_grid(plot_rw_a,legend,rel_widths = c(3,.25))
      }
      if(par=='b'){ ###rw beta=====
          if(form=='stan'){
            post=rstan::extract(mod)
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
          theme_classic(14)+
          theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                strip.text = element_text(face="bold", size=12),
                axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
        
        
        legend = cowplot::get_legend(plot1)
        
        plot_rw_b=cowplot::plot_grid(plot1 + theme(legend.position="none"),
                                     plot2 + theme(legend.position="none"),
                                     ncol=2,nrow=1,labels=c("A","B"))
        plot=cowplot::plot_grid(plot_rw_b,legend,rel_widths = c(3,.25))
      }
      if(par=='both'){ #rw alpha beta=====
        if(form=='stan'){
          post=rstan::extract(mod)
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
            geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2]),linewidth=1.3)+
            geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3]),linewidth=1.3)+
            geom_point(aes(colour = gamma_df$gamma),size=2.5) +
            scale_colour_viridis_c(name='p')+
            ggtitle(title)+
            xlab("Spawners") + 
            ylab("Recruits")+
            theme_classic(14)+
            theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                  strip.text = element_text(face="bold", size=12),
                  axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
          
          
          
          legend = cowplot::get_legend(plot1)
          
          plot_hmm_a=cowplot::plot_grid(plot1 + theme(legend.position="none"),
                                       plot2 + theme(legend.position="none"),
                                       ncol=2,nrow=1,labels=c("A","B"))
          plot=cowplot::plot_grid(plot_hmm_a,legend,rel_widths = c(3,.3))
          
        }
        
        if(par=='b'){ #hmm beta====
          
          if(form=='stan'){
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
            geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2]),linewidth=1.3)+
            geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3]),linewidth=1.3)+
            geom_point(aes(colour = gamma_df$gamma),size=2.5) +
            scale_colour_viridis_c(name='p')+
            ggtitle(title)+
            xlab("Spawners") + 
            ylab("Recruits")+
            theme_classic(14)+
            theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                  strip.text = element_text(face="bold", size=12),
                  axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
          
         
          
          legend = cowplot::get_legend(plot1)
          
          plot_hmm_b=cowplot::plot_grid(plot1 + theme(legend.position="none"),
                                        plot2 + theme(legend.position="none"),
                                        ncol=2,nrow=1,labels=c("A","B"))
          plot=cowplot::plot_grid(plot_hmm_b,legend,rel_widths = c(3,.3))
        }
        if(par=='both'){ #hmm alpha beta====
          if(form=='stan'){
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
            geom_line(data=pred_df,aes(x=x_new,y=pred_df[,2]),linewidth=1.3)+
            geom_line(data=pred_df,aes(x=x_new,y=pred_df[,3]),linewidth=1.3)+
            geom_point(aes(colour = gamma_df$gamma),size=2.5) +
            scale_colour_viridis_c(name='p')+
            ggtitle(title)+
            xlab("Spawners") + 
            ylab("Recruits")+
            theme_classic(14)+
            theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
                  strip.text = element_text(face="bold", size=12),
                  axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
          
          legend = cowplot::get_legend(plot1)
          
          plot_hmm_ab=cowplot::plot_grid(plot1 + theme(legend.position="none"),
                                        plot2 + theme(legend.position="none"),
                                        ncol=2,nrow=1,labels=c("A","B"))
          plot=cowplot::plot_grid(plot_hmm_ab,legend,rel_widths = c(3,.3))
        }
    }
    if(make.pdf==TRUE){
      pdf(here(path,paste(paste(title,type,par,form,sep='_'),'.pdf',sep='')))
      return(plot)
      dev.off()
    }
  if(make.pdf==FALSE){
    return(plot)
  }
  
  }
  