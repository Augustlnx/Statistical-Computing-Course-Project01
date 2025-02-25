set.seed(42)
realpa=0.6
realpb=0.4
E=function(nk,a,b,alpha){
  A=(choose(5,nk)*a^nk*(1-a)^(5-nk))
  B=(choose(5,nk)*b^nk*(1-b)^(5-nk))
  Evalue=alpha*A/(alpha*A+(1-alpha)*B)
  return(Evalue)
}
K=100
max=200

#
#finalbootpa=c()
#finalbootpb=c()
#finaljackpa=c()
#finaljackpb=c()
finalbootpa_mse=c()
finalbootpb_mse=c()
finaljackpa_mse=c()
finaljackpb_mse=c()
#

for(z in 1:7){

samplesize=z*40-30
B=500
bootpa=c()
bootpb=c()
jackpa=c()
jackpb=c()
bootpa_mse=c()
bootpb_mse=c()
jackpa_mse=c()
jackpb_mse=c()

for(i in 1:K){
  #产生样本
  n=c()
  for(j in 1:samplesize){
    prob <- ifelse(runif(1) < 0.5, realpa, realpb)
    samples <-rbinom(5, 1, prob)
    n[j] <- sum(samples)
  }
  #BOOTSTRAP
  bstpa = numeric(B)
  bstpb = numeric(B)
  bstpa_mse = numeric(B)
  bstpb_mse = numeric(B)
  for(j in 1:B){
    m=sample(seq(1,samplesize),samplesize,replace=T)
    bootn=n[m]
    pa = numeric(max + 1)
    pb = numeric(max + 1)
    alpha = numeric(max + 1)
    pa[1]=0.6
    pb[1]=0.4
    alpha[1]=0.1
    for(l in 1:max){
      EY=E(bootn,pa[l],pb[l],alpha[l])
      pa[l+1]=sum(bootn*EY)/(5*sum(EY))
      pb[l+1]=sum(bootn*(1-EY))/(5*sum(1-EY))
      alpha[l+1]=sum(EY)/samplesize
      if(abs(pa[l+1]-pa[l])<1e-8&abs(pb[l+1]-pb[l])<1e-8&abs(alpha[l+1]-alpha[l])<1e-8) break
    }
    bstpa[j]<-pa[length(pa)]
    bstpb[j]<-pb[length(pb)]
    #MSE
    bstpa_mse[j]<-(bstpa[j]-realpa)^2
    bstpb_mse[j]<-(bstpb[j]-realpb)^2
  }
  #bootpa[i]<-mean(bstpa)
  #bootpb[i]<-mean(bstpb)
  bootpa_mse[i]<-mean(bstpa_mse)
  bootpb_mse[i]<-mean(bstpb_mse)
  #JACKKNIFE
  jkpa = numeric(samplesize)
  jkpb = numeric(samplesize)
  jkpa_mse=numeric(samplesize)
  jkpb_mse=numeric(samplesize)
  for(j in 1:samplesize){
    jackn=n[-j]
    pa = numeric(max + 1)
    pb = numeric(max + 1)
    alpha = numeric(max + 1)
    pa[1]=0.6
    pb[1]=0.4
    alpha[1]=0.1
    for(l in 1:max){
      EY=E(jackn,pa[l],pb[l],alpha[l])
      pa[l+1]=sum(jackn*EY)/(5*sum(EY))
      pb[l+1]=sum(jackn*(1-EY))/(5*sum(1-EY))
      alpha[l+1]=sum(EY)/(samplesize-1)
      if(abs(pa[l+1]-pa[l])<1e-8&abs(pb[l+1]-pb[l])<1e-8&abs(alpha[l+1]-alpha[l])<1e-8) break
    }
    jkpa[j]<-pa[length(pa)]
    jkpb[j]<-pb[length(pb)]
    #MSE
    jkpa_mse[j]<-(jkpa[j]-realpa)^2
    jkpb_mse[j]<-(jkpb[j]-realpb)^2
  }
  #jackpa[i]<-mean(jkpa)
  #jackpb[i]<-mean(jkpb)
  jackpa_mse[i]<-mean(jkpa_mse)
  jackpb_mse[i]<-mean(jkpb_mse)
}

#finalbootpa[z]=mean(bootpa)
#finalbootpb[z]=mean(bootpb)
#finaljackpa[z]=mean(jackpa)
#finaljackpb[z]=mean(jackpb)
finalbootpa_mse[z]=mean(bootpa_mse)
finalbootpb_mse[z]=mean(bootpb_mse)
finaljackpa_mse[z]=mean(jackpa_mse)
finaljackpb_mse[z]=mean(jackpb_mse)

}

dataPA <- data.frame(
  accuracy = c(bootpa, jackpa),
  method = factor(c(rep("Bootstrap", length(bootpa)), rep("Jackknife", length(jackpa))))
)

# 使用ggplot2绘制密度图，并添加红色垂直线
p_plot <- ggplot(dataPA, aes(x = accuracy, fill = method)) +
  geom_density(alpha = 0.5, adjust = 1) + # adjust参数用于控制曲线的平滑度
  geom_vline(xintercept = 0.6, color = "red", linetype = "solid", size = 1) + # 添加红色垂直线
  labs(
    title = "Density of PA Estimates for Bootstrap and Jackknife",
    x = "P",
    y = "Density",
    fill = "Method" # 添加图例标签
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("#42BCB2", "#235689")) + # 手动设置填充颜色
  theme(legend.title = element_blank() )

ggsave(filename = paste0("density_plot_PA_z", z, ".png"), plot = p_plot, width = 8, height = 6)

dataPB <- data.frame(
  accuracy = c(bootpb, jackpb),
  method = factor(c(rep("Bootstrap", length(bootpb)), rep("Jackknife", length(jackpb))))
)

# 使用ggplot2绘制密度图，并添加红色垂直线
pb_plot <- ggplot(dataPB, aes(x = accuracy, fill = method)) +
  geom_density(alpha = 0.5, adjust = 1) + # adjust参数用于控制曲线的平滑度
  geom_vline(xintercept = 0.4, color = "red", linetype = "solid", size = 1) + # 添加红色垂直线
  labs(
    title = "Density of PB Estimates for Bootstrap and Jackknife",
    x = "P",
    y = "Density",
    fill = "Method" # 添加图例标签
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("#42BCB2", "#235689")) + # 手动设置填充颜色
  theme(legend.title = element_blank())

ggsave(filename = paste0("density_plot_PB_z", z, ".png"), plot = pb_plot, width = 8, height = 6)

  

#随样本量增加的变化图
#pa
n = c(10, 50, 90, 130, 170, 210,250)
datanpa <- data.frame(
  n = c(10, 50, 90, 130, 170, 210,250),
  bootstrap_estimate = finalbootpa,
  jackknife_estimate = finaljackpa
)

# 绘制图表
ggplot(datanpa, aes(x = n)) +
  geom_line(aes(y = bootstrap_estimate, color = "Bootstrap"), size = 1) +
  geom_line(aes(y = jackknife_estimate, color = "Jackknife"), size = 1) +
  geom_point(aes(y = bootstrap_estimate, color = "Bootstrap"), size = 3) +
  geom_point(aes(y = jackknife_estimate, color = "Jackknife"), size = 3) +
  geom_hline(yintercept = 0.6, color = "black", linetype = "dashed", size = 1) +
  # 使用annotate()来添加文本标注
  annotate("text", x = max(datanpa$n) + 10, y = 0.6, label = "True value: 0.6", hjust = 0) +
  labs(title = "Estimation of PA with Increasing Sample Size (B=500) ",
       x = "Sample Size (n)",
       y = "Estimated PA",
       color = "Method") +
  theme_minimal() +
  theme(legend.title = element_blank())

#pb
datanpb <- data.frame(
  n = c(10, 50, 90, 130, 170, 210,250),
  bootstrap_estimate = finalbootpb,
  jackknife_estimate = finaljackpb
)

# 绘制图表
ggplot(datanpb, aes(x = n)) +
  geom_line(aes(y = bootstrap_estimate, color = "Bootstrap"), size = 1) +
  geom_line(aes(y = jackknife_estimate, color = "Jackknife"), size = 1) +
  geom_point(aes(y = bootstrap_estimate, color = "Bootstrap"), size = 3) +
  geom_point(aes(y = jackknife_estimate, color = "Jackknife"), size = 3) +
  geom_hline(yintercept = 0.4, color = "black", linetype = "dashed", size = 1) +
  # 使用annotate()来添加文本标注
  annotate("text", x = max(datanpb$n) + 10, y = 0.4, label = "True value: 0.4", hjust = 0) +
  labs(title = "Estimation of PB with Increasing Sample Size (B=500) ",
       x = "Sample Size (n)",
       y = "Estimated PB",
       color = "Method") +
  theme_minimal() +
  theme(legend.title = element_blank())

install.packages("htmlwidgets")









#MSE柱状图
n_values <- c(10, 50, 90, 130, 170, 210, 250)
bootpa_mse <- finalbootpa_mse 
jackpa_mse <- finaljackpa_mse


# 创建数据框
dat <- data.frame(
  n = rep(n_values, times = 2),  # 每个样本量重复2次
  method = rep(c("Bootstrap PA", "Jackknife PA"), each = length(n_values)),  # 定义方法
  mse = c(finalbootpa_mse, finaljackpa_mse)  # 合并两种方法的MSE值
)

# 绘制分组条形图，并添加美化元素
plot_ly(dat,
        x = ~n, color = ~method, y = ~mse,
        type = "bar", barmode = "group",
        marker = list(line = list(color = "black", width = 1.5)),
        text = ~paste0("样本量：", n, "<br>", "MSE：", round(mse, 2)),
        hoverinfo = "text"
) %>%
  layout(
    title = "Bootstrap PA与Jackknife PA的MSE比较（B=500）",
    xaxis = list(
      title = "样本量",
      tickvals = n_values,
      ticktext = paste0(n_values, "个样本")
    ),
    yaxis = list(
      title = "MSE",
      tickformat = ".2f"
    ),
    legend = list(
      title = "方法",
      x = 0.85,  # 调整图例位置
      y = 0.95,
      bgcolor = "rgba(255, 255, 255, 0.8)",  # 图例背景色
      bordercolor = "black",
      borderwidth = 1
    ),
    bargap = 0.15,  # 条形之间的间隔
    bargroupgap = 0.1,  # 组之间的间隔
    plot_bgcolor = "white",  # 图表背景色
    paper_bgcolor = "white",  # 页面背景色
    margin = list(l = 50, r = 50, t = 70, b = 50)  # 图表边距
  ) %>%
  config(displayModeBar = FALSE)

