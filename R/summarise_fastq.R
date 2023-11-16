#' A function to summarise the fastq from a run.
#'
#' @description
#' This function is heavy and depends a lot on the depth of your fastq. Be careful when setting the cores, this function can saturate the RAM.
#'
#'
#'
#' @param fastq Object returned by \link{list_fastq} function.
#' @param cores The number of cores to be used. Default = 10.
#' @param number_of_fastq Number of fastq to integrate in the summary. Default = 10
#'
#' @return
#' @export
#'
#' @examples
#' # See Vignette.
summarise_fastq<- function(fastq, cores= 1, number_of_fastq=10){
require(ShortRead)
require(parallel)
theme_set(theme_classic()+
  theme(plot.title = element_text(size = cores, face="bold"),
  axis.title = element_text(size=20, face="bold"),
  axis.text = element_text(size=20)))


fastq= rng_fastq(fastq, n=number_of_fastq)

if(length(fastq)>2){
  print("pair end analysis")
fwd <- mclapply(X= fastq$fastq_fwd, readFastq, mc.cores = cores)
rev = mclapply(X= fastq$fastq_rv, readFastq, mc.cores = cores)

width_fwd <- purrr::map(fwd, ShortRead::width)%>%
  purrr::map(sample, size=1000)%>%
  unlist

width_rev <- purrr::map(rev, ShortRead::width)%>%
  purrr::map(sample, size=1000)%>%
  unlist

df_fwd= purrr::map(fwd, function(x){
  tmp= as(quality(x), 'matrix')
  data.frame(median= apply(tmp, 2, median, na.rm=T),
             mean = apply(tmp, 2, mean, na.rm=T),
             variance= apply(tmp, 2, function(x)sd(x, na.rm=T)/mean(x, na.rm=T)),
             bp= 1:dim(tmp)[2])
  })

names(df_fwd)= fastq$names

df_rev= purrr::map(rev, function(x){
  tmp= as(quality(x), 'matrix')
  data.frame(median= apply(tmp, 2, median, na.rm=T),
             mean = apply(tmp, 2, mean, na.rm=T),
             variance= apply(tmp, 2, function(x)sd(x, na.rm=T)/mean(x, na.rm=T)),
             bp= 1:dim(tmp)[2])
})

names(df_rev)= fastq$names

df_fwd= do.call("rbind", df_fwd)
df_rev= do.call("rbind", df_rev)

p_fwd= df_fwd%>%
  group_by(bp)%>%
  mutate(median=mean(median),
         variance=mean(variance),
         mean=mean(mean))%>%
  ggplot()+
  geom_line(aes(y=median/100, bp, color="Median"), linewidth=1)+
  geom_line(aes(y=variance, bp, color="Variance"), linetype=2, linewidth=1)+
  geom_point(aes(y=mean/100, bp, color="Mean"), size=2)+
  geom_density(data= as.data.frame(width_fwd), aes(width_fwd, after_stat(count/sum(count)), fill="Density"), adjust=.10, alpha=0.3, position="identity")+
  scale_colour_manual("", breaks = c("Median", "Variance", "Mean"), values = c("orange", "darkgreen", "salmon"))+
  scale_fill_manual("", breaks = c("Density"), values=c("grey"))+
scale_y_continuous(labels = function(x)x*100)+
  labs(x="Base",
y="Quality Score")
p_rev= df_rev %>%
  group_by(bp) %>%
  mutate(median=mean(median),
         variance=mean(variance),
         mean=mean(mean))%>%
  ggplot()+
  geom_line(aes(y=median/100, bp, color="Median"), linewidth=1)+
  geom_line(aes(y=variance, bp, color="Variance"), linetype=2, linewidth=1)+
  geom_point(aes(y=mean/100, bp, color="Mean"), size=2)+
  geom_density(data= as.data.frame(width_fwd), aes(width_fwd, after_stat(count/sum(count)), fill="Density"), adjust=.10, alpha=0.3, position="identity")+
  scale_colour_manual("", breaks = c("Median", "Variance", "Mean"), values = c("orange", "darkgreen", "salmon"))+
  scale_fill_manual("", breaks = c("Density"), values=c("grey"))+
  scale_y_continuous(labels = function(x)x*100)+
  labs(x="Base",
       y="Quality Score")
p= ggpubr::ggarrange(p_fwd, p_rev, nrow = 1, common.legend = T)

df = list(depth= list(fwd=as.numeric( summary(fwd)[,1]), rev =as.numeric(summary(rev)[,1])), # length of each sample
          quality= list(fwd= df_fwd, rev=df_rev),
          length=list(fwd= sample(width_fwd, 1000), rev= sample(width_rev, 1000)),
          sample_name= c(paste(fastq$names, "fwd", sep="_"), paste(fastq$names, "rev", sep="_")))

  } else {

  print("single end analysis")
fwd <- mclapply(X= fastq$fastq_fwd, readFastq, mc.cores = cores)
width_fwd <- purrr::map(fwd, ShortRead::width)%>%
  # purrr::map(sample, size=1000)%>%
  unlist

df_fwd= purrr::map(fwd, function(x){
  tmp= as(quality(x), 'matrix')
  data.frame(median= apply(tmp, 2, median, na.rm=T),
             mean = apply(tmp, 2, mean, na.rm=T),
             variance= apply(tmp, 2, function(x)sd(x, na.rm=T)/mean(x, na.rm=T)),
             bp= 1:dim(tmp)[2]
)
})
names(df_fwd)= fastq$names

df_fwd= do.call("rbind", df_fwd)

p= df_fwd%>%
  group_by(bp)%>%
  mutate(median=mean(median),
         variance=mean(variance),
         mean=mean(mean))%>%
  ggplot()+
  geom_line(aes(y=median/100, bp, color="Median"), linewidth=1)+
  geom_line(aes(y=variance, bp, color="Variance"), linetype=2, linewidth=1)+
  geom_point(aes(y=mean/100, bp, color="Mean"), size=2)+
  geom_density(data= as.data.frame(width_fwd), aes(width_fwd, after_stat(count/sum(count)), fill="Density"), adjust=.10, alpha=0.3, position="identity")+
  scale_colour_manual("", breaks = c("Median", "Variance", "Mean"), values = c("orange", "darkgreen", "salmon"))+
  scale_fill_manual("", breaks = c("Density"), values=c("grey"))+
  scale_y_continuous(labels = function(x)x*100)+
  labs(x="Base", y="Quality Score")

df = list(depth= list(fwd=as.numeric( summary(fwd)[,1]), rev =as.numeric(summary(fwd)[,1])), # length of each sample
quality= list(fwd= df_fwd),
length=list(fwd= summary(width_fwd)["Mean"]),
sample_name= c(paste(fastq$names, "fwd", sep="_")))



  }
return(list(plot= p, data= df))
}
