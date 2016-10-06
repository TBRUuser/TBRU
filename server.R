library(shiny)
library(ggvis)
library(dplyr)
library(ggplot2)
library(cowplot)
library(plotly)
library(DT)
library(reshape2)


#####################################--------------| 
##### Dataset Input Requirement #####              |
#####################################              |
## Id                                              |
## dilution_or (original data of dilution)         |
## Sample (peptide)                                |
## IFNg                                            |
## ------------------------------------------------|


#######################
##### Preparation #####
#######################
cv <- function(x) {
  mm = mean(x)
  if(mm == 0) return(0)
  res = sqrt(var(x)) / mm
  res
}


###########################
##### Data Processing #####
###########################


## Heatmap ##
breaks1=c(0,500,1000)                       ## the scale of value, (low,mid,high)
colkey1=c("red","yellow","green")           ## corresponding color


########################
##### Shiny Server #####
########################

shinyServer(
  function(input,output,session){
    
    ##################################################################
    ## 1.1 Introduction ## --- had been replaced by other function in UI.R
    #    output$Intro1<-renderPlot(
    #      {
    #        try<-readImage("C:/Users/Horace/OneDrive/Project/TB/T1/picture/introduction.jpg")
    #        display(try,method = "raster")
    #      }
    #    )
    
    
    ######################################################
    ########## Data Entry & Download ##########
    ######################################################
    data<-reactive({
      if (is.null(input$dataset)) {return()}
      switch(input$dataset,
             "dataset0.csv" = read.csv("./data/dataset0.csv"),
             "dataset2.csv" = read.csv("./data/dataset2.csv"),
             "dataset3.csv" = read.csv("./data/dataset3.csv")
      )
    })
    
    output$download_data<-downloadHandler(
      filename = function() {paste0("./data/",input$dataset)},
      content = function(file) {write.csv(data(),file)}
    )
    
    ##################################################################
    ###################################
    ##### ggvis tooltips function #####
    ###################################
    ### x in the function means the whole row (an individual in dataset) ###
    
    ## use in 1.2.2.4  DE - plot - line plot ##
    vis_tooltip <- function(x) {
      if (is.null(x)) return(NULL)
      paste0(
        names(x),":",format(x),collapse = "<br />"
      )
    }
    
    
    ##################################################################
    ## 1.2.1.1  DE - Summary - Dataset, Show Raw Data ##
    output$rawdata <- DT::renderDataTable({
      data<-data()
      datatable(data)
    })
    
    
    ##################################################################
    ##  1.2.1.2 -- 1.2.1.3  DE - Summary - table & Distribution ##
    output$tb_vis_dilu<-renderTable(
      {
        data<-data()
#        p<-table(data$Visit, data$Dilution)       ## original version, must work
        p<-data.frame(data$Visit, data$Dilution)
        names(p)=c("Visit","Dilution")
        table(p)
      }
    )
    output$tb_sample<-renderTable(
      {
        data<-data()
#        table(data$Sample)    ## original version, must work.
        Sample<-data.frame(data$Sample)   ## this line have risk, be careful
        names(Sample)=c("Sample")
        table(Sample)
      }
    )
    output$histo1<-renderPlotly(
      {
        data<-data()
        mean.well = aggregate(IFNg ~ ID + Dilution+Sample+Visit, data, mean)
        cv.well = aggregate(IFNg ~ ID + Dilution+Sample+Visit, data, cv)
        histo1<-plot_ly(data=cv.well,x=cv.well$IFNg,type = "histogram",opacity=0.8)
        histo1<-layout(histo1,xaxis=list(title="cross validation of IFNg"),yaxis=list(title="Frequency"))
#        hist(cv.well$IFNg,10,main="Histogram of IFNg",xlab="IFNg after adjustment by cv.well",ylab="Frequency")
      }
    )
    output$histo2<-renderPlotly(
      {
        data<-data()
        mean.well = aggregate(IFNg ~ ID + Dilution+Sample+Visit, data, mean)
        histo2<-plot_ly(data=mean.well,x=mean.well$IFNg,type = "histogram",opacity=0.8)
        histo2<-layout(histo2,xaxis=list(title="IFNg"),yaxis=list(title="Frequency"))
      }
    )
    
    ##################################################################
    ## 1.2.1.4  DE - Summary - CV analysis ##
    output$box1<-renderPlotly(
      {
        data<-data()
        ## compute variability among wells ##
        mean.well = aggregate(IFNg ~ ID + Dilution+Sample+Visit, data, mean)
        cv.well = aggregate(IFNg ~ ID + Dilution+Sample+Visit, data, cv)
        ii = which( data$ID==311 & data$Dilution==2 & data$Sample=="UNS" & data$Visit==1)
        ##  inter-individual variability
        mm0 = aggregate(IFNg ~ ID + Dilution+Sample+Visit, data, mean)
        cv.ind = aggregate(IFNg ~ Dilution+Sample+Visit, mm0, cv)
        mean.ind = aggregate(IFNg ~ Dilution+Sample+Visit, mm0, mean)
        ## inter-visit variability
        cv.visit = aggregate(IFNg ~ ID + Dilution+Sample, mm0, cv)
        mean.visit = aggregate(IFNg ~ ID + Dilution+Sample, mm0, mean)
        allcv1 = list(cv.well$IFNg, cv.ind$IFNg, cv.visit$IFNg)
        names(allcv1) = c("Well", "Individual", "Visit")
        ## exclude items with mean 0
        ix=mean.well[,5]>5
        cv.well.no0 = cv.well[ix,5]*100
        ix=mean.ind[,4]>5
        cv.ind.no0 = cv.ind[ix,4]*100
        ix=mean.visit[,4]>5
        cv.visit.no0 = cv.visit[ix,4]*100
        allcv2 = list(cv.well.no0, cv.ind.no0, cv.visit.no0)
        names(allcv2) = c("Well", "Individual", "Visit")
        allcv3 = list(cv.visit.no0, cv.ind.no0)
        names(allcv3) = c("Assay", "Individual")
        ## from here down I'll collapse wells
        data0 = mm0
        ## reshape list table for interactive plot_ly ##
        allcv1_m<-melt(allcv1)
        allcv2_m<-melt(allcv2)
        allcv3_m<-melt(allcv3)
        ## interactive plot axis ##
        xa<-list(title="")
        ya<-list(title="IFNg")
        ## interactive boxplot ##
        if (input$box=="inter-1")
        {
          inter1<-plot_ly(data=allcv1_m,x = allcv1_m$L1,y = allcv1_m$value,color = allcv1_m$L1,
                          type = "box")
          inter1<-layout(inter1,xaxis=xa,yaxis=ya,
                         title="Cross Validation of IFNg (among well, ID or Visit)")   ## well = all
        }
        else if (input$box=="inter-2") 
        {
          inter2<-plot_ly(data=allcv2_m,x = allcv2_m$L1,y = allcv2_m$value,color = allcv2_m$L1,
                          type = "box")
          inter2<-layout(inter2,xaxis=xa,yaxis=ya,
                         title="Cross Validation of IFNg (among well, ID or Visit), Exclude Item with Mean 0")   ## well = all
        }
        else if (input$box=="inter-3")
        {
          inter3<-plot_ly(data=allcv3_m,x = allcv3_m$L1,y = allcv3_m$value,color = allcv3_m$L1,type = "box")
          inter3<-layout(inter3,xaxis=xa,yaxis=ya)
        }
       else NULL
        
#        if (input$box=="inter-1")
#        {boxplot(allcv1, ylab="Coefficient of variation of IFNg")}
#        else if (input$box=="inter-2") 
#        {boxplot(allcv2, ylab="Coef. var. (mean/sd), in %")}
#        else if (input$box=="inter-3")
#        {boxplot(allcv3, ylab="Variability (%)")}
#        else NULL
      }
    )
    
    ##################################################################
    ## 1.2.2.1  DE - plot - Boxplot ##
    ############# Old Ver Boxplot -- fast but not interactive #############    
    #    output$box2<-renderPlot(
    #      {
    #        box2ex<-data[data$Dilution==input$dilution_b2 & data$Visit==input$visit_b2,]
    #        plot(x=box2ex$Sample,y=box2ex$IFNg,las=2,main="Boxplot of Peptides' IFNg with Fixed Dilution & Visit",ylab="IFNg")
    #      }
    #    )
    
    output$ui_dilution_box2<-renderUI({
      data<-data()
      selectInput("dilution_b2","Dilution value:",levels(factor(data$Dilution)))
    })
    
    output$ui_visit_box2<-renderUI({
      data<-data()
      selectInput("visit_b2","Visit:",levels(factor(data$Visit)))
    })
    
    output$box2<-renderPlotly({
      data<-data()
      box2ex<-data[data$Dilution==input$dilution_b2 & data$Visit==input$visit_b2,]
      vbox2<-plot_ly(data=box2ex,x=box2ex$Sample,y=box2ex$IFNg,type="box")
      margin<-list(l=80,b=200)
      xa_vbox2<-list(title="Sample")
      ya_vbox2<-list(title="IFNg")
      layout(vbox2,margin=margin,xaxis=xa_vbox2,yaxis=ya_vbox2)
    })
    
    
    ## 1.2.2.2  DE - plot - Variety across Visit ##
    
    #    output$plot3<-renderPlot(
    #      {
    #        plot3ex<-data[data$Dilution==input$dilution_vav & data$ID==input$ID_vav,]
    #        plot(x=plot3ex$Sample,y=plot3ex$IFNg,las=2,main="Boxplot of Peptides' IFNg with Fixed Dilution & ID",ylab="IFNg")
    #      }
    #    )
    
    output$ui_dilution_vav<-renderUI({
      data<-data()
      selectInput("dilution_vav","Dilution value:",levels(factor(data$Dilution)))
    })
    
    output$ui_ID_vav<-renderUI({
      data<-data()
      selectInput("ID_vav","Present ID  (across PEPTIDE):",levels(factor(data$ID)))
    })
    
    output$plot3<-renderPlotly({
      data<-data()
      plot3ex<-data[data$Dilution==input$dilution_vav & data$ID==input$ID_vav,]
      bplot3<-plot_ly(x=plot3ex$Sample,y=plot3ex$IFNg,type="box",
                      main="Boxplot of Peptides' IFNg with Fixed Dilution & ID",ylab="IFNg")
      margin<-list(l=80,b=200)
      xa_bplot3<-list(title="Sample")
      ya_bplot3<-list(title="IFNg")
      layout(bplot3,margin=margin,xaxis=xa_bplot3,yaxis=ya_bplot3)
    })
    
    
    ##################################################################
    ## 1.2.2.3  DE - plot - Heatmap ##
    
    ### First one - Plotly interavtice plot ###
    output$ui_dilution_hmp1<-renderUI({
      data<-data()
      selectInput("dilution_hmp1","Dilution value:",levels(factor(data$Dilution)))
    })
    
    output$ui_visit_hmp1<-renderUI({
      data<-data()
      selectInput("visit_hmp1","Visit:",levels(factor(data$Visit)))   
    })
    
    output$heatmp1<-renderPlotly(
      {
        data<-data()
        mm0 = aggregate(IFNg ~ ID + Dilution+Sample+Visit, data, mean)
        data0 = mm0
        ##   breaks1=c(0,500,1000)  &  colkey1=c("red","yellow","green")   adjust in   "box_heat.R"   ##
        heatmp1ex<-data0[data0$Dilution==input$dilution_hmp1 & data0$Visit==input$visit_hmp1,]
        heatmp1ex<-data.frame(heatmp1ex$ID,heatmp1ex$Sample,heatmp1ex$IFNg)
        colnames(heatmp1ex)<-c("ID","Sample","IFNg")
        heatmp1ex$ID<-factor(heatmp1ex$ID)
        gg1<-ggplot(heatmp1ex,aes(y=ID,x=Sample,fill=IFNg));
        gg1<-gg1 + geom_tile(color="white", size=0.1);
        gg1<-gg1 + scale_fill_gradientn(name="IFNg",colors = colkey1,breaks=breaks1,labels=format(breaks1));
        gg1<-gg1 + geom_text(aes(fill=heatmp1ex$IFNg,label=round(heatmp1ex$IFNg,0)),size=3);
        gg1<-gg1 + xlab("Sample") + ylab("ID");
        gg1<-gg1 + theme(axis.text.x=element_text(angle = 70,vjust = 1,hjust = 1));
        gg1<-ggplotly(gg1);
        margin<-list(l=80,b=200)
        layout(gg1,margin=margin)
      }
    )
    
    ### Second One - ggvis - update from usual ggplot ###
    output$ui_dilution_hmp2<-renderUI({
      data<-data()
      selectInput("dilution_hmp2","Dilution value:",levels(factor(data$Dilution)))
    })
    
    output$ui_visit_hmp2<-renderUI({
      data<-data()
      selectInput("visit_hmp2","Visit:",levels(factor(data$Visit)))   
    })
    
    output$heatmp2<-renderPlot(
      {
        data<-data()
        mm0 = aggregate(IFNg ~ ID + Dilution+Sample+Visit, data, mean)
        data0 = mm0
        ##   breaks1=c(0,500,1000)  &  colkey1=c("red","yellow","green")   adjust in   "box_heat.R"   ##
        heatmp2ex<-data0[data0$Dilution==input$dilution_hmp2 & data0$Visit==input$visit_hmp2,]
        heatmp2ex<-data.frame(heatmp2ex$ID,heatmp2ex$Sample,heatmp2ex$IFNg)
        colnames(heatmp2ex)<-c("ID","Sample","IFNg")
        heatmp2ex$ID<-factor(heatmp2ex$ID)
        gg2<-ggplot(heatmp2ex,aes(y=ID,x=Sample,fill=IFNg))
        gg2<-gg2 + geom_tile(color="white", size=0.1)
        gg2<-gg2 + scale_fill_gradientn(name="IFNg",colors = colkey1,breaks=breaks1,labels=format(breaks1))
        gg2<-gg2 + geom_text(aes(fill=heatmp2ex$IFNg,label=round(heatmp2ex$IFNg,0)))
        gg2<-gg2 + xlab("Sample") + ylab("ID")
        gg2<-gg2 + theme(axis.text.x=element_text(angle = 70,vjust = 1,hjust = 1)) 
        #gg2<-gg2 + coord_equal()
        gg2<-ggdraw(switch_axis_position(gg2,axis = "x"))
        print(gg2)
      }
    )
    
    
    ##################################################################
    ## 1.2.2.4  DE - plot - line plot - overall ##
    output$ui_title_lpall<-renderUI({
      if (is.null(input$across_type)) {return()}
      switch (input$across_type,
        "Peptides" = h3("Overall Peptides' IFNg change among visit"),
        "Individuals" = h3("Overall ID's IFNg among visit")
      )
    })
    
    output$ui_dilution_lpall<-renderUI({
      if (is.null(input$across_type)) {return()}
      data<-data()
      switch (input$across_type,
        "Peptides" = selectInput("dilution_lpap","Dilution value:",levels(factor(data$Dilution))),
        "Individuals" = selectInput("dilution_lpai","Dilution value:",levels(factor(data$Dilution)))
      )
    })
    
    output$ui_others_lpall<-renderUI({
      if (is.null(input$across_type)) {return()}
      data<-data()
      switch (input$across_type,
        "Peptides" = selectInput("ID_lpap","Present ID  (across PEPTIDE):",levels(factor(data$ID))),
        "Individuals" = selectInput("peptides_lpai","Peptites  (across INDIVIDUAL):",levels(factor(data$Sample)))
      )
    })
    
    output$ui_plot_lpall<-renderUI({
      if (is.null(input$across_type)) {return()}
      data<-data()
      switch (input$across_type,
        "Peptides" = ggvisOutput("plot2"),
        "Individuals" = ggvisOutput("plot1")
      )
    })
    
    ## 1.2.2.4  DE - plot - line plot - across individual - Overall ##
    vis1<-reactive({
      data<-data()
      data0 = aggregate(IFNg ~ ID + Dilution+Sample+Visit, data, mean)
      dat.dilution = data0[data0$Dilution==input$dilution_lpai,]
      thisData = dat.dilution[dat.dilution$Sample == input$peptides_lpai, c(1,4,5)]
      thisData %>% ggvis(~Visit,~IFNg,stroke=~factor(ID)) %>%
        layer_points() %>%
        layer_lines() %>%
        add_axis("x",value=c(levels(factor(data$Visit)))) %>%
        add_legend("stroke",title="ID") %>%
#        hide_legend("stroke") %>%
        add_tooltip(vis_tooltip,"hover") 
    })
    vis1 %>% bind_shiny("plot1")
    
    ## 1.2.2.4  DE - plot - line plot - across peptide - Overall ##
    vis2<-reactive({
      data<-data()
      data0 = aggregate(IFNg ~ ID + Dilution+Sample+Visit, data, mean)
      dat.dilution = data0[data0$Dilution==input$dilution_lpap,]
      thisData = dat.dilution[dat.dilution$ID == input$ID_lpap, c(3,4,5)]
      thisData %>% ggvis(~Visit,~IFNg,stroke=~factor(Sample)) %>%
        layer_points()  %>%
        layer_lines()  %>%
        add_axis("x",value=c(levels(factor(data$Visit)))) %>%
        add_legend("stroke",title="Sample") %>%
#        hide_legend("stroke") %>%
        add_tooltip(vis_tooltip,"hover") 
    })
    vis2 %>% bind_shiny("plot2")
    
    
    ##################################################################
    ## 1.2.2.4  DE - plot - line plot - single ##
    output$ui_title_lpsingle<-renderUI({
      if (is.null(input$across_type)) {return()}
      switch (input$across_type,
        "Peptides" = h3("Single Peptides' IFNg change among visit"),
        "Individuals" = h3("Single ID's IFNg among visit")
      )
    })
    
    output$ui_dilution_lpsingle<-renderUI({
      if (is.null(input$across_type)) {return()}
      data<-data()
      switch (input$across_type,
              "Peptides" = selectInput("dilution_lpap_s","Dilution value:",levels(factor(data$Dilution))),
              "Individuals" = selectInput("dilution_lpai_s","Dilution value:",levels(factor(data$Dilution)))
      )
    })
    
    output$ui_others_lpsingle<-renderUI({
      if (is.null(input$across_type)) {return()}
      data<-data()
      switch (input$across_type,
              "Peptides" = selectInput("ID_lpap_s","Present ID  (across PEPTIDE):",levels(factor(data$ID))),
              "Individuals" = selectInput("peptides_lpai_s","Peptites  (across INDIVIDUAL):",levels(factor(data$Sample)))
      )
    })
    
    output$ui_checkbox_lpsingle<-renderUI({
      if (is.null(input$across_type)) {return()}
      data<-data()
      switch (input$across_type,
        "Peptides" = checkboxGroupInput("peptides_lpap_s_check","Peptides",c(levels(factor(data$Sample)))),
        "Individuals" = checkboxGroupInput("ID_lpai_s_check","ID",c(levels(factor(data$ID))))
      )
    })
    
    output$ui_plot_lpsingle<-renderUI({
      if (is.null(input$across_type)) {return()}
      data<-data()
      switch (input$across_type,
              "Peptides" = plotlyOutput("lpap_s",height = "500px"),
              "Individuals" = plotlyOutput("lpai_s",height = "500px")
      )
    })
    
    #######################
    ## Select All Button ##
    observe({
      data<-data()
      if (input$across_type=="Peptides"){
        updateCheckboxGroupInput(
          session,"peptides_lpap_s_check",choices=c(levels(factor(data$Sample))),
          selected = if (input$select_all) c(levels(factor(data$Sample)))
       )
      } else if (input$across_type=="Individuals") {
        updateCheckboxGroupInput(
          session,"ID_lpai_s_check",choices=c(levels(factor(data$ID))),
          selected = if (input$select_all) c(levels(factor(data$ID)))
       )
      }
    })
    
    ###########################################
    ## Base Ver by ggplot2 - across peptides ##
    output$lpap_s<-renderPlotly({
      data<-data()
      data0 = aggregate(IFNg ~ ID + Dilution+Sample+Visit, data, mean)
      dat.dilution = data0[data0$Dilution==input$dilution_lpap_s,]
      thisdata_all<-dat.dilution
      thisdata_all<-filter(thisdata_all,ID==input$ID_lpap_s)
      thisdata<-NULL
      if (is.null(input$peptides_lpap_s_check)) {return()}
      else {
        for (i in 1:length(input$peptides_lpap_s_check)) {
          thisdata_piece<-filter(thisdata_all,Sample==input$peptides_lpap_s_check[i])
          thisdata<-rbind(thisdata,thisdata_piece)
        }
        ### can change the ggplot to ggvis from here, for interactive ###
        gg3<-ggplot(data=thisdata,aes(x=Visit,y=IFNg,colour=Sample))+geom_point()
        if (is.numeric(thisdata$Visit[1])) {
#          gg3<-gg3+geom_line(aes(group=thisdata$Sample))+scale_x_continuous(breaks = as.numeric(levels(factor(thisdata$Visit))))
          gg3<-gg3+geom_line()+scale_x_continuous(breaks = as.numeric(levels(factor(thisdata$Visit))))
        }
        else {
#          gg3<-gg3+geom_line(aes(group=thisdata$Sample))
          gg3<-gg3+geom_line()
        }
        gg3<-gg3+xlab("Visit")+ylab("IFNg")
        gg3<-gg3+scale_color_discrete("Sample")
        gg3<-gg3+theme(plot.title=element_text(size=24))
        gg3<-ggplotly(gg3)
        margin<-list(l=80,r=100,b=100)
        layout(gg3,margin=margin)  
      }
    })
    
    ##############################################
    ## Base Ver by ggplot2 - across individuals ##
    output$lpai_s<-renderPlotly({
      data<-data()
      data0 = aggregate(IFNg ~ ID + Dilution+Sample+Visit, data, mean)
      dat.dilution = data0[data0$Dilution==input$dilution_lpai_s,]
      thisdata_all<-dat.dilution
      thisdata_all<-filter(thisdata_all,Sample==input$peptides_lpai_s)
      thisdata<-NULL
      if (is.null(input$ID_lpai_s_check)) {ggplot()}
      else {
        for (i in 1:length(input$ID_lpai_s_check)) {
          thisdata_piece<-filter(thisdata_all,ID==input$ID_lpai_s_check[i])
          thisdata<-rbind(thisdata,thisdata_piece)
        }
        ### can change the ggplot to ggvis from here, for interactive ###
        gg4<-ggplot(data=thisdata,aes(x=Visit,y=IFNg,colour=factor(ID)))+geom_point()
        if (is.numeric(thisdata$Visit[1])) {
#          gg4<-gg4+geom_line(aes(group=factor(thisdata$ID)))+scale_x_continuous(breaks = as.numeric(levels(factor(thisdata$Visit))))
          gg4<-gg4+geom_line()+scale_x_continuous(breaks = as.numeric(levels(factor(thisdata$Visit))))
        }
        else {
          gg4<-gg4+geom_line()
        }
        gg4<-gg4+xlab("Visit")+ylab("IFNg")
        gg4<-gg4+scale_color_discrete("ID")
        gg4<-gg4+theme(plot.title=element_text(size=24))
        gg4<-ggplotly(gg4)
        margin<-list(l=80,r=100,b=100)
        layout(gg4,margin=margin)  
      }
    })
    
    ## Base Ver by ggvis - across peptides ##  
    #    cbs1<-reactive({
    #      dat.dilution = data0[data0$Dilution==input$dilution_tc,]
    #      thisdata_all<-dat.dilution
    #      thisdata_all<-filter(thisdata_all,ID==input$ID_tc)
    #      thisdata<-NULL
    #      if (is.null(input$line_tc)) {thisdata<-thisdata_all}
    #      else {
    #        for (i in 1:length(input$line_tc)) {
    #          thisdata_piece<-filter(thisdata_all,Sample==input$line_tc[i])
    #          thisdata<-rbind(thisdata,thisdata_piece)
    #        }
    #      }
    #      thisdata %>% ggvis(~Visit,~IFNg,stroke=~factor(Sample)) %>%
    #        layer_points() %>%
    #        layer_lines() %>%
    #        add_axis("x",value=c(levels(factor(data$Visit)))) %>%
    #        add_legend("stroke",title="Sample") %>%
    #        hide_legend("stroke") %>%
    #        add_tooltip(vis_tooltip,"hover") 
    #    })
    #    cbs1 %>% bind_shiny("check_test_plot")

    
    ##################################################################
    
    
  }
)