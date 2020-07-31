{ 
  library(ropls)
  library(ggplot2)
  library(ggrepel)
  library(FactoMineR)
  library(ComplexHeatmap)
  library(readxl)
  library("circlize")
  library(DT)
  library(shiny)
  library(shinybootstrap2) 
  library(corrplot)
  library(plotly)
  library(car)
}
IQRdata<-function(x) {
  
  
  Q_range=1;Q_shang=1;Q_xia=1#这个是弄个初始值，否则下面回报错
  colx=ncol(x)
  rowx=nrow(x)
  
  while(sum(boxplot(t(x),plot=F)[["out"]])!=0)
  {
    for (i in 1:rowx) {
      q<-x[i,]
      q<-as.matrix(q[q!=0])
      Q_xia[[i]] <- quantile(q, probs = 0.25)
      Q_shang[[i]] <- quantile(q, probs = 0.75)
      Q_range[[i]] <- Q_shang[[i]] - Q_xia[[i]]
      for (j in 1:colx) {
        if(!is.na(x[i,j])){
        if (x[i,j] > (Q_shang[[i]] + 1.5*Q_range[[i]]))
        {
          x[i,j]<-0
        }
        if (x[i,j] < (Q_xia[[i]] - 1.5*Q_range[[i]]))
        {
          x[i,j]<-0
        }
      }}}}
  return(x)
}
s<-function(x,group)
{x.x = t(x)
k.x = matrix(group, ncol = 1)
x.n = cbind(k.x, x.x)
sorted = x.n[order(x.n[, 1]), ]
k = matrix(sorted[, 1], ncol = 1)
g = c()
for (i in 1:nrow(sorted)) {
  if (any(g == sorted[i, 1])) {
    g = g
  }
  else {
    g = matrix(c(g, sorted[i, 1]), ncol = 1)
  }
}
Y = matrix(rep(NA, nrow(sorted)), ncol = 1)
for (i in 1:nrow(sorted)) {
  for (l in 1:2) {
    if (sorted[i, 1] == l) {
      Y[i, ] = 0
    }
    else {
      Y[i, ] = 1
    }
  }
}
X = as.matrix(sorted[, -1], ncol = ncol(sorted) - 1)
nf = 1
T = c()
P = c()
C = c()
W = c()
Tortho = c()
Portho = c()
Wortho = c()
Cortho = c()
for (j in 1:nf) {
  w = (t(X) %*% Y) %*% solve(t(Y) %*% Y)
  w1 = t(w) %*% w
  w2 = abs(sqrt(w1))
  w = w %*% solve(w2)
  t = (X %*% w) %*% solve(t(w) %*% w)
  t1 = t(t) %*% t
  c = t(Y) %*% t %*% solve(t1)
  c1 = t(c) %*% c
  u = Y %*% c %*% solve(c1)
  u1 = t(u) %*% u
  u2 = abs(sqrt(u1))
  p = (t(X) %*% t) %*% solve(t1)
  wortho = p - w
  wortho1 = t(wortho) %*% wortho
  wortho2 = abs(sqrt(abs(wortho1)))
  wortho = wortho %*% solve(wortho2)
  tortho = X %*% wortho %*% solve(t(wortho) %*% wortho)
  tortho1 = t(tortho) %*% tortho
  portho = t(X) %*% tortho %*% solve(tortho1)
  cortho = t(Y) %*% tortho %*% solve(tortho1)
  X = X - tortho %*% t(portho)
  T = matrix(c(T, t))
  P = matrix(c(P, p))
  C = matrix(c(C, c))
  W = matrix(c(W, w))
  Tortho = matrix(c(Tortho, tortho))
  Portho = matrix(c(Portho, portho))
  Wortho = matrix(c(Wortho, wortho))
  Cortho = matrix(c(Cortho, cortho))
}
T = matrix(T, ncol = nf)
T = scale(T, scale = FALSE, center = TRUE)

s = as.matrix(sorted[, -1], ncol = ncol(sorted) - 1)
p1 = c()
for (i in 1:ncol(s)) {
  scov = cov(s[, i], T)
  p1 = matrix(c(p1, scov), ncol = 1)
}
pcorr1 = c()
for (i in 1:nrow(p1)) {
  den = apply(T, 2, sd) * sd(s[, i])
  corr1 = p1[i, ]/den
  pcorr1 = matrix(c(pcorr1, corr1), ncol = 1)
}
p1<<-p1
pcorr1<<-pcorr1
ss<-data.frame(p1,pcorr1)
return(ss)
}

paretoscale<-function(z){
  rowmean<-apply(z,1,mean)
  rowsd<-apply(z,1,sd)
  rowsqrtsd<-sqrt(rowsd)
  rv<-sweep(z,1,rowmean,"-")
  rv<-sweep(rv,1,rowsqrtsd,"/")
  return(rv)
}
ui <-{ fluidPage(
  sidebarPanel(
    fileInput('file1', '点击一下Browse',
              accept=c('text/csv', 'text/comma-separated-values,text/plain')),
    tags$hr(),
    checkboxInput('header', 'Header', TRUE),
    actionButton("submit", "确定怎么改"),
    uiOutput("fj"),
    uiOutput("tl") ,
    uiOutput("wk") ,
    downloadLink('downloadData', 'Download'),
    
    downloadButton("Data", "分析结果
                   (vip,p,FC和标准化的数据)"),
    sliderInput(inputId = "fz",
                label = "分多少组",
                min=0,
                max=20,
                value = 2
    ),
    selectInput(inputId = "del",
                label = "删除异常值",
                choices = c("IQR","NO"),
                selected = "NO"),
    sliderInput(inputId = "del2",
                label = "删除缺失值为多少的行",
                min=0,
                max=20,
                value = 5
    ),
    selectInput(inputId = "bu",
                label = "补缺失值",
                choices = c("最小值的1/2","每一列的平均值","NO"),
                selected = "最小值的1/2"),
    selectInput(inputId = "bzh",
                label = "标准化",
                choices = c("总离子流标准化+log化","平均值标准化","NO"),
                selected = "总离子流标准化+log化"),
    h1("下面进行分析"),
    
    selectInput(inputId = "pc",
                label = "pca标准化",
                choices = c("uv格式化","ctr格式化","Par标准化","NO"),
                selected = "uv格式化 "),
    selectInput(inputId = "pls",
                label = "(o)pls标准化",
                choices = c("uv格式化","ctr标准化","Par标准化","standard标准化","NO"),
                selected = "uv格式化"),
    sliderInput(inputId = "plsda",
                label = "(o)pls正交次数(0次即plsda)",
                min=0,
                max=10,
                value = 1
    ),
    numericInput(inputId = "hst",
                 label = "火山图log2FC阈值",
                 value=1),
    numericInput(inputId = "hst2",
                 label = "火山图p阈值",
                 value=0.05),
    numericInput(inputId = "vip",
                 label = "热图vip阈值",
                 value=1)
  ),
  mainPanel(
    DT::dataTableOutput('contents',height = 40),
    
   mainPanel(
      tabsetPanel(type = "tabs",
                  #tabPanel()语句可以使主面板像侧边栏面板一样工作
                  tabPanel("3DPCA", br(),plotlyOutput("pca", width = 800, height = 800)),
                  tabPanel("PCA", br(),plotOutput("mpgPlot", width = 800, height = 800)),
                  tabPanel("OPLS", br(),plotOutput("opls", width = 800, height = 800)),
                  tabPanel("Splot", br(), plotOutput("splot", width = 800, height = 800)),
                  tabPanel("火山图", br(),plotOutput("hst", width = 800, height = 800)),
                  tabPanel("热图", br(), plotOutput("rt", width = 800, height = 800)),
                  tabPanel("相关性", br(), plotOutput("xg", width = 700, height = 800))
                 ),
      
    ),
   
   )
    
  )
}

server <- function(input, output) {
  url1 <- a("富集分析", href="https://www.metaboanalyst.ca/MetaboAnalyst/upload/PathUploadView.xhtml") 
  output$fj <- renderUI({ 
    tagList("把你的差异代谢物上传到:", url1) 
  })
  url2 <- a("kegg通路分析", href="https://www.kegg.jp/kegg/tool/map_pathway2.html") 
  output$tl <- renderUI({ 
    tagList("然后把keggID上传到:", url2) 
  })
  url3 <- a("悟空平台", href="https://www.kegg.jp/kegg/tool/map_pathway2.html") 
  output$wk <- renderUI({ 
    tagList("或者直接上传到这个平台:", url3) 
  })
  qcl <- reactive({
    {inFile <- input$file1
    ctr<-(read.csv(inFile$datapath, header=input$header,row.names = 1))
    name<-ctr[,1]
    ctr<-ctr[,-1]
    if(input$del=="IQR"){ ctr<-IQRdata(ctr)}
    if(input$bu=="最小值的1/2"){v=min(ctr[ctr!=0])/2;ctr[ctr==0]<-v}
    ctr<-as.data.frame(ctr)
    if(input$bzh=="总离子流标准化+log化"){
      ctr<-sweep(ctr,2,apply(ctr, 2, sum),FUN = "/")
    }
    }
    ctr<-data.frame(name,ctr)
    return(ctr)
  })
  output$contents <-DT::renderDataTable((ctr<-(read.csv(input$file1$datapath, header=input$header))), filter = 'top', options = list( pageLength = 5, autoWidth = TRUE))
  qclopls <- reactive({
    ctr<-as.data.frame(qcl())
    name<-ctr[,1]
    ctr<-ctr[,-1]
    if(input$bzh=="总离子流标准化+log化"){
      ctr<-log(ctr,2)
    }
    ctr<-as.matrix(ctr)
    
    d<-"none"
    if(input$pls=="uv格式化"){ctr<-arm::rescale(ctr,"UV")}
    if(input$pls=="ctr格式化"){d<-"center"}
    if(input$pls=="Par标准化"){d<-"pareto"}
    if(input$pls=="standard标准化"){d<-"standard"}
    
    ctr<-as.data.frame(ctr)
    group<-rep(c(paste("C",1:input$fz)),each=(ncol(ctr)/input$fz))##分组的
    aa <-opls(t(ctr), group,predI = 1, orthoI = input$plsda,crossvalI =ncol(ctr) ,
              scaleC=d)
    aa<<-aa
    return(aa)
  })
  qclp <- reactive({
    ctr<-qcl()[,-1]
    if(input$bzh=="总离子流标准化+log化"){
      ctr<-log(ctr+1,2)
    }
    b<-nrow(ctr)
    ctr<-as.data.frame(ctr)
    padj<-1
    f<-(ncol(ctr)/input$fz)
    for (i in 1:b) {
      padj[i]<-t.test(c(as.numeric(ctr[i,1:f])),c(as.numeric(ctr[i,(1+f):(2*f)])))["p.value"]
    }
    padj<-(padj)
    return(padj)
  })
  qclFC <- reactive({
    ctr<-qcl()[,-1]
    ctr<<-as.data.frame(ctr)
    f<-(ncol(ctr)/input$fz)
    F1<-t(ctr[,1:f]);F2<-t(ctr[,(1+f):(2*f)])
    FC<-colMeans(F1)/ colMeans(F2)
    return(FC)
  })
  output$hst <- renderPlot({
    ctr<-(qcl())
    FC<-qclFC()
    log2FC<-log(FC,2)
    padj<-as.numeric(as.data.frame(qclp()))
    data <- data.frame(log2FC=log2FC,padj=padj,qcl()[,1])#和并在一起
    colnames(data)[3]<-"name"
    data$col<-"no"#设置分类，先全部假设为没有no
    data$col[data$padj <= input$hst2 & data$log2FC >= input$hst] <- "up"#如果padj值小于0.05以及log2FC大于等于2，就让数据的col那一列对应的数据no变成up
    data$col[data$padj <= input$hst2 & data$log2FC <= -input$hst] <- "down"#如果padj值小于0.05以及log2FC小于等于-2，就让数据的col那一列对应的数据no变成down
    data<-as.data.frame(data)
    xlim=max(-log2FC,log2FC) #求出xlim的绝对值的最大值，其实是xlim的最大值和-xlim的最大值比较，输出两者最那个
    #下面将画图的过程存在p中
    p<-ggplot(data=data,aes(x=log2FC,y=-1*log10(padj),color = col))+#我有一数据集名字叫data,它的x坐标的数据为log2FC，y坐标数据-1*log10(padj),颜色用我的col数据集
      geom_point()+ #请将它用散点图画出来
      labs(x="log2(FC)",y="-log10(FDR)")+ #它x轴名称为log2(FC)，y轴名称为-log10(FDR)
      ggtitle("Volcano plot")+#它的标题名称为Volcano plot
      scale_color_manual(values =c("#619cff" ,"grey" ,"#f8766d"))+
      geom_hline(aes(yintercept=-1*log10(input$hst2)),colour="black", linetype="dashed") +#加一条y=-1*log10(0.05)的直线，它的颜色为black，线条类型为dashed即虚线
      geom_vline(xintercept=c(-input$hst,input$hst),colour="black", linetype="dashed")+#加一条x=2，和一条x=-2的直线，它的颜色为black，线条类型为dashed即虚线
      theme_bw()+#将背景设置为theme_bw类型，然后再储存在volcano中
      geom_label_repel(data = data[data$col!="no",],#将图的的data数据集的col不等于no和vip大于1的选出来
                       aes(label = name))#将它的名字标上去
    print(p)
  })
  output$pca <- renderPlotly({
    name<-qcl()[,1]
    ctr<-qcl()[,-1]
    ctr<-as.matrix(ctr)
    if(input$pc=="uv格式化"){ctr<-arm::rescale(ctr,"UV")}
    if(input$pc=="ctr格式化"){ctr<-arm::rescale(ctr,"center")}
    if(input$pc=="Par标准化"){ctr<-paretoscale(ctr)}
    ctr<-as.data.frame(ctr)
    group<-rep(c(paste("C",1:input$fz)),each=(ncol(ctr)/input$fz))##分组的
    pca<-PCA(t(ctr),graph = FALSE)
    dat_1 <- as.data.frame(pca$ind$coord) #dat_1为主成分得分（样本坐标）
    name<-rownames(dat_1)
    dat_1<-cbind(dat_1,name)
    pc1_ctr <- paste(round(pca$eig[1,2], 2), "%", sep = "") #PC1的方差贡献率
    pc2_ctr <- paste(round(pca$eig[2,2], 2), "%", sep = "") #PC2的方差贡献率
    pc3_ctr <- paste(round(pca$eig[3,2], 2), "%", sep = "") #PC2的方差贡献率
    
    
    fig <- plot_ly(dat_1, x = ~Dim.1, y = ~Dim.2, z = ~Dim.3, color = group, colors = c('#BF382A', '#0C4B8E'))
    fig <- fig %>% add_markers()
    fig <- fig %>% layout(scene = list(xaxis = list(title = pc1_ctr),
                                       yaxis = list(title = pc2_ctr),
                                       zaxis = list(title = pc3_ctr)))
  })
  output$mpgPlot <- renderPlot({
    name<-qcl()[,1]
    ctr<-qcl()[,-1]
    ctr<-as.matrix(ctr)
    if(input$pc=="uv格式化"){ctr<-arm::rescale(ctr,"UV")}
    if(input$pc=="ctr格式化"){ctr<-arm::rescale(ctr,"center")}
    if(input$pc=="Par标准化"){ctr<-paretoscale(ctr)}
    ctr<-as.data.frame(ctr)
    group<-rep(c(paste("C",1:input$fz)),each=(ncol(ctr)/input$fz))##分组的
    pca<-PCA(t(ctr),graph = FALSE)
    dat_1 <- as.data.frame(pca$ind$coord) #dat_1为主成分得分（样本坐标）
    name<-rownames(dat_1)
    dat_1<-cbind(dat_1,name)
    pc1_ctr <- paste(round(pca$eig[1,2], 2), "%", sep = "") #PC1的方差贡献率
    pc2_ctr <- paste(round(pca$eig[2,2], 2), "%", sep = "") #PC2的方差贡献率
    #dat_2 <- as.data.frame(pca$var$coord) #dat_2为主成分载荷（变量坐标）
    ggplot(dat_1, aes(Dim.1, Dim.2,color = group)) + #映射PC1,PC2
      geom_point(size = 7) + #绘制得分散点图
      labs(title = "PCA") + #图表标题
      xlab(paste("PC1", pc1_ctr, sep = "  ")) + #x轴标题
      ylab(paste("PC2", pc2_ctr, sep = "  ")) + #y轴标题
      stat_ellipse( geom = "polygon", alpha = 0.1, level = 0.95, show.legend = F,
                    col="black",size=1)+#加置信椭圆
      geom_vline(xintercept=mean(dat_1$Dim.1),lty=1,col="black",lwd=0.5) +#添加横线
      geom_hline(yintercept =mean(dat_1$Dim.2),lty=1,col="black",lwd=0.5)+#添加竖线
      theme( plot.background = element_rect(fill="white"),
             panel.background = element_rect(fill='white', colour='gray'),
             strip.text.x=element_text(size=rel(1.2), family="serif", angle=-90),
             strip.text.y=element_text(size=rel(1.2),  family="serif") ,
             axis.text.x = element_text(size = 14,color="black"),
             axis.text.y = element_text(size = 20,color="black"),axis.ticks.x=element_blank()
      )+
      geom_text_repel(data = dat_1 ,
                      aes(label = name),
                      size = 5,
                      check_overlap = TRUE,
                      colour = "red",
                      box.padding = unit(0.2, "lines"),
                      point.padding = unit(0.1, "lines"), segment.color = "black", show.legend = FALSE )
  })
  output$opls <- renderPlot({
    plot(qclopls()) 
  }) 
  output$splot<- renderPlot({
    ctr<-as.data.frame(qcl())
    name<-ctr[,1]
    ctr<-ctr[,-1]
    if(input$bzh=="总离子流标准化+log化"){
      ctr<-log(ctr,2)
    }
    ctr<-as.matrix(ctr)
    if(input$pls=="uv格式化"){ctr<-arm::rescale(ctr,"UV")}
    if(input$pls=="ctr格式化"){ctr<-arm::rescale(ctr,"center")}
    if(input$pls=="Par标准化"){ctr<-paretoscale(ctr)}

    ctr<-as.data.frame(ctr)
    group<-rep(c(1:input$fz),each=(ncol(ctr)/input$fz))##分组的
    plot(s(ctr,group))
  }) 
  da <- reactive({
      dat<-data.frame(qclopls()@vipVn,as.numeric(as.data.frame(qclp())) ,qclFC(),qcl())
       colnames(dat)<-c("VIP","p","FC",colnames(qcl()))
       x<-as.data.frame(dat)
       dat<<-x
    return(x)
  })
  output$Data <- downloadHandler(
    filename = function() {
      paste('分析结果-', Sys.Date(), '.csv', sep='')
    },
    content = function(file) {
      dat<-data.frame(qclopls()@vipVn,as.numeric(as.data.frame(qclp())) ,qclFC(),qcl())
      colnames(dat)<-c("VIP","p","FC",colnames(qcl()))
      write.csv(dat, file)
    }
  )
  output$downloadData <- downloadHandler(
     filename = function() {
        paste('分析结果-', Sys.Date(), '.csv', sep='')
      },
      content = function(con) {
        write.csv(da() , con)
      }
    )

  output$rt <- renderPlot({
    dat<-data.frame(qclopls()@vipVn,as.numeric(as.data.frame(qclp())) ,qclFC(),qcl())
    colnames(dat)<-c("VIP","p","FC","name",colnames(qcl()[,-1]))
    x<-as.data.frame(dat)
    x$FC<-log(x$FC,2)
    x<-x[((x$VIP>input$vip)&(x$FC>input$hst|x$FC<(-(input$hst)))&(x$p<input$hst2)&!is.na(x$name)),]
    row.names(x)<-x$name
    x<-x[,5:ncol(x)]
    group<-rep(c(paste("C",1:input$fz)),each=(ncol(qcl()[,-1])/input$fz))##分组的
    f<-(ncol(ctr)/input$fz)
    FC<-colMeans(t(x[,1:f]));
    FC1<-colMeans(t(x[,(1+f):(2*f)]))
    FC<-data.frame(FC,FC1)
    FC<<-as.matrix(FC)
    col=colorRamp2(c(-2,0,2),c("green","white","red"))
    x<-t(scale(t(x)))
    x<<-as.matrix(x)
    ha = rowAnnotation(foo = anno_lines(FC,  width = unit(2, "cm"),
                                        gp = gpar(col = c(5,4)), 
                                        add_points = TRUE, 
                                        pt_gp = gpar(col = c(5,4)),
                                        pch = c(1,16)))
    d<-Heatmap(x,
               top_annotation = HeatmapAnnotation(group=anno_block(gp = gpar(fill = c(5,4)),
                                                                   labels = c("C", "H"), 
                                                                   labels_gp = gpar(col = "white", fontsize = 10))),
               left_annotation = ha,col =col,column_split = factor(rep(c(paste("C",1:input$fz)),each=(ncol(qcl()[,-1])/input$fz))),
               column_names_gp = gpar(col = rep(c(5,4)),each=(ncol(qcl()[,-1])/input$fz),fontsize=20))
    
    print(d)
  output$xg <- renderPlot({
    dat<-data.frame(qclopls()@vipVn,as.numeric(as.data.frame(qclp())) ,qclFC(),qcl())
    colnames(dat)<-c("VIP","p","FC",colnames(qcl()))
    x<-as.data.frame(dat)
    x$FC<-log(x$FC,2)
    x<-x[(x$VIP>input$vip)&(x$FC>input$hst|x$FC<(-(input$hst)))&(x$p<input$hst2)&!is.na(x$name),]
    rownames(x)<-x$name
    x<-x[,5:ncol(x)]
    M <- cor(t(x))
    corrplot(M, method = "circle")
    
    }) 
  })
}
shinybootstrap2::withBootstrap2({shinyApp(ui, server)})
#library(miniUI)
#runGadget(ui, server) #打开应用
bv<-data.frame(dat[,c(1,4)],p1,pcorr1)
write.csv(bv,"data.csv")
plot(p1,pcorr1)
x<-read.csv("data.csv")
plot(x$p1,x$pcorr1)
x$col<-"red"
x[x$p1>0.1,]
ggplot(data=x)+geom_point(aes(p1,x$pcorr1,col=col))+
geom_label_repel(data = x[x$p>0,],#将图的的data数据集的col不等于no和vip大于1的选出来
                 aes(label = name))#将它的名字标上去
x[order(x$p1),][(nrow(x)-9):nrow(x),]<-"red"
