pdf(file='mouse_output.pdf',width=4.5,height=4.5);
gstable=read.table('mouse_output.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("H3c1","H2bc11","H2ac10","Btbd35f1","Gm20738","H3c13","Mfap1a","Gm14150","Gm10058","Gm10439")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='Rep1,Rep2,Rep3_vs_plasmid neg.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(474.94034284360606,46.196621525103915,17.50023984017576,48.23002373913933),c(771.9452896359455,561.7839153320673,739.6768039114287,600.8316516655493),c(477.61606308497846,25.573129772825382,10.500143904105455,4.904748176861627),c(551.1983697227203,179.83684807986882,212.33624339413254,203.5470493397575),c(271.5856044993015,8.249396700911413,5.833413280058586,9.809496353723254),c(559.2255304468375,79.19420832874957,56.000767488562424,51.49985585704708),c(651.5378787741863,177.3620290695954,145.83533200146465,227.2533321945887),c(339.8164706542984,122.09107117348891,63.000863424632726,203.5470493397575),c(387.97943499900214,8.249396700911413,10.500143904105455,19.61899270744651),c(389.3172951196883,23.098310762551957,66.50091139266787,40.05544344436995),c(369.2493933093951,13.199034721458261,17.50023984017576,22.07136679587732),c(291.6535063095947,47.021561195195055,60.667498112609294,49.864939798093204),c(678.2950811879106,122.09107117348891,74.6676899847499,95.64258944880173),c(370.5872534300813,419.8942920763909,361.6716233636323,457.77649650708514),c(472.2646226022336,86.61866535956985,36.16716233636323,96.46004747827865),c(414.7366374127264,136.11504556503832,143.50196668944122,140.6027810700333),c(280.950625344105,7.424457030820272,2.3333653120234343,12.261870442154066),c(381.29013439557104,36.29734548401022,11.666826560117173,13.896786501107943),c(548.5226494813478,39.597104164374784,9.333461248093737,36.7856113264622),c(551.1983697227203,41.246983504557065,26.833701088269496,45.77764965070852))
targetgene="H3c1"
collabel=c("plasmid","Rep1","Rep2","Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(659.5650394983036,14.848914061640544,35.00047968035152,25.34119891378507),c(603.3749144294826,29.697828123281088,17.50023984017576,21.253908766400382),c(284.9642057061636,86.61866535956985,116.66826560117171,165.94397998381837),c(680.970801429283,25.573129772825382,12.833509216128888,19.61899270744651),c(351.8572117404743,14.848914061640544,21.00028780821091,12.261870442154066),c(548.5226494813478,11.549155381275979,2.3333653120234343,2.4523740884308136),c(413.39877729204017,13.199034721458261,19.83360515219919,10.626954383200191),c(465.5753219988025,51.146259545650764,40.833892960410104,48.23002373913933),c(516.4140065848786,73.41963063811158,134.16850544134746,58.85697812233952))
targetgene="H2bc11"
collabel=c("plasmid","Rep1","Rep2","Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(595.3477537053653,238.40756465633984,504.0069073970618,236.24537051883502),c(533.8061881537994,34.647466143827934,87.50119920087879,40.87290147384689),c(495.0082446538992,41.246983504557065,14.000191872140606,18.80153467796957),c(476.27820296429223,588.1819847749838,288.1706160348941,315.53879937809796),c(313.0592682405741,317.60177298508944,112.00153497712485,189.65026283864958),c(369.2493933093951,34.647466143827934,39.66721030439838,49.047481768616265),c(417.4123576540988,28.872888453189947,21.00028780821091,27.793573002215886),c(504.3732654987027,405.87031768484155,400.172151012019,479.8478633029625),c(240.81482172351855,158.38841665749914,292.837346658941,238.69774460726583),c(426.77737849890235,87.44360502966099,56.000767488562424,66.21410038763196),c(519.0897268262511,25.573129772825382,29.167066400292928,19.61899270744651),c(374.60083379213995,38.77216449428364,49.00067155249212,64.57918432867808),c(159.2053543616595,108.89203645203065,70.00095936070304,66.21410038763196),c(513.7382863435063,80.84408766893185,39.66721030439838,89.10292521298622),c(426.77737849890235,108.06709678193951,93.33461248093738,92.37275733089398),c(509.7247059814476,377.8223689017427,353.5048447715503,470.85582497871616),c(440.1559797057645,60.220595916653316,25.667018432257777,47.412565709662395),c(517.7518667055649,47.846500865286195,30.333749056304647,46.59510768018546),c(505.71112561938895,70.11987195774701,64.16754608064444,68.66647447606277),c(382.6279945162573,136.11504556503832,178.50244636979272,297.5547227296054),c(529.7926077917408,61.8704752568356,102.66807372903111,85.01563506560153))
targetgene="H2ac10"
collabel=c("plasmid","Rep1","Rep2","Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(346.50577125772946,61.04553558674446,74.6676899847499,106.26954383200192),c(417.4123576540988,147.6642009463143,53.66740217653899,194.55501101551118),c(573.9419917743859,454.54175822021887,497.0068114609915,635.9823469330576),c(604.7127745501688,156.73853731731685,201.83609949002707,179.02330845544938),c(687.6601020327141,340.7000837476414,238.00326182639031,326.9832117907751),c(587.3205929812481,414.94465405584407,350.00479680351515,344.96728843926775),c(722.4444651705556,329.1509283663654,319.6710477472105,289.380142434836),c(559.2255304468375,104.76733810157495,126.00172684926545,154.49956757114126),c(418.750217774785,333.2756267168211,267.17032822668324,153.68210954166432),c(317.07284860263275,273.0550308001678,513.3403686451555,434.070213652254),c(181.94897641332514,11.549155381275979,29.167066400292928,118.53141427415598),c(738.4987866187902,424.8439300969378,558.8409922296125,675.2203323479506),c(442.8316999471369,438.86790448848717,415.3390255401713,433.25275562277704),c(290.3156461889085,272.23009113007663,189.0025902738982,199.45975919237281),c(181.94897641332514,29.697828123281088,51.334036864515554,44.96019162123158),c(369.2493933093951,123.7409505136712,93.33461248093738,73.5712226529244),c(309.0456878785155,142.71456292576744,268.3370108826949,398.1020603552687),c(501.6975452573303,174.88721005932197,316.17099977917536,246.05486687255828),c(405.3716165679229,85.79372568947869,162.1688891856287,139.78532304055636),c(398.6823159644918,290.3787638720818,64.16754608064444,210.08671357557301),c(565.9148310502686,196.33564148169162,259.0035496346012,319.62608952548266),c(473.60248272291983,207.88479686296762,149.3353799694998,166.7614380132953),c(404.03375644723667,354.7240581391908,313.83763446715193,455.3241224186543),c(581.9691524985031,219.43395224424359,85.16783388885536,217.44383584086546),c(659.5650394983036,248.30684069743353,341.8380182114331,302.459470906467),c(551.1983697227203,146.83926127622317,221.66970464222626,218.2612938703424),c(412.060917171354,390.19646395310986,283.50388541084726,157.76939968904898),c(398.6823159644918,127.04070919403577,189.0025902738982,218.2612938703424),c(519.0897268262511,413.2947747156618,536.6740217653899,350.6894946456063))
targetgene="Btbd35f1"
collabel=c("plasmid","Rep1","Rep2","Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(418.750217774785,53.621078555924186,74.6676899847499,104.63462777304804),c(604.7127745501688,429.79356811748465,686.0094017348897,391.5623961194532),c(311.7214081198879,381.1221275821073,421.1724388202299,419.3559691216691),c(362.56009270596405,185.6114257705068,129.5017748173006,255.0469051968046),c(393.33087548174694,17.323733071913967,214.66960870615597,96.46004747827865),c(409.38519692998153,183.13660676023338,135.33518809735918,298.3721807590823),c(509.7247059814476,140.23974391549402,368.6717192997026,295.9198066706515),c(490.9946642918406,537.035725229333,833.0114163923661,538.704841425302),c(347.8436313784157,324.20129034581856,395.5054203879721,351.50695267508326),c(408.0473368092953,221.908771254517,277.6704721307887,167.57889604277224),c(339.8164706542984,156.73853731731685,137.66855340938264,122.61870442154067),c(341.1543307749846,140.23974391549402,46.66730624046869,71.93630659397053),c(611.4020751535999,445.4674218492163,390.83868976392523,297.5547227296054),c(252.85556280969448,70.94481162783815,38.500527648386665,96.46004747827865),c(628.7942567225207,348.9494804485528,281.17052009882383,307.36421908332863),c(276.93704498204636,187.26130511068908,58.334132800585856,120.16633033310985),c(307.70782775782925,113.01673480248637,85.16783388885536,65.39664235815502),c(438.81811958507825,57.745776906379895,73.50100732873818,128.34091062787923),c(539.1576286365444,111.36685546230407,219.33633933020283,83.38071900664765),c(421.4259380161575,320.901531665454,575.1745494137766,406.2766406500381),c(366.5736730680227,235.93274564606642,285.8372507228707,253.41198913785072),c(464.23746187811633,99.81770008102811,166.83561980967556,118.53141427415598),c(856.230477239177,660.7766757430043,659.1757006466202,724.2678141165669))
targetgene="Gm20738"
collabel=c("plasmid","Rep1","Rep2","Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(489.65680417115436,202.93515884242078,77.00105529677333,96.46004747827865),c(476.27820296429223,14.848914061640544,24.50033577624606,25.34119891378507),c(331.7893099301811,30.52276779337223,51.334036864515554,16.349160589538755),c(32.10864289646914,5.774577690637989,4.666730624046869,2.4523740884308136),c(489.65680417115436,273.0550308001678,194.83600355395677,200.27721722184975),c(448.18314042988175,256.55623739834493,106.16812169706627,112.80920806781742),c(572.6041316536997,124.56589018376233,113.16821763313656,103.8171697435711),c(283.6263455854774,44.546742184921634,3.5000479680351515,21.253908766400382),c(42.811523861958854,42.896862844739346,9.333461248093737,13.896786501107943),c(319.7485688440052,20.623491752278532,44.33394092844525,26.15865694326201),c(446.8452803091955,118.79131249312435,60.667498112609294,108.72191792043273))
targetgene="H3c13"
collabel=c("plasmid","Rep1","Rep2","Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(232.78766099940125,12.37409505136712,22.166970464222626,13.079328471631005),c(307.70782775782925,10.724215711184836,18.666922496187475,14.71424453058488),c(326.43786944743624,42.071923174648205,40.833892960410104,45.77764965070852),c(421.4259380161575,68.46999261756473,127.16840950527717,69.48393250553971),c(325.10000932675007,58.570716576471035,25.667018432257777,37.60306935593914),c(373.2629736714538,327.5010490261831,80.50110326480848,58.85697812233952),c(418.750217774785,144.36444226594972,106.16812169706627,93.19021536037091))
targetgene="Mfap1a"
collabel=c("plasmid","Rep1","Rep2","Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(490.9946642918406,86.61866535956985,100.33470841700768,99.72987959618641),c(656.8893192569311,35.472405813919075,54.83408483255071,36.7856113264622),c(437.480259464392,28.872888453189947,18.666922496187475,13.896786501107943),c(433.4666791023334,15.673853731731686,11.666826560117173,30.2459470906467),c(404.03375644723667,67.64505294747359,229.83648323430828,56.40460403390871),c(492.3325244125268,137.7649249052206,52.50071952052727,82.56326097717071))
targetgene="Gm14150"
collabel=c("plasmid","Rep1","Rep2","Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(529.7926077917408,30.52276779337223,79.33442060879676,30.2459470906467),c(337.14075041292597,190.56106379105364,110.83485232111313,146.32498727637187),c(710.4037240843797,405.0453780147504,522.6738298932493,506.82397827570145),c(335.80289029223974,44.546742184921634,51.334036864515554,138.1504069816025),c(299.680667033712,16.498793401822827,128.33509216128888,23.706282854831198),c(325.10000932675007,108.89203645203065,71.16764201671475,127.5234525984023),c(406.70947668860913,102.29251909130153,120.16831356920687,192.9200949565573),c(382.6279945162573,194.68576214150934,126.00172684926545,105.45208580252498),c(545.8469292399753,216.95913323397016,235.66989651436685,168.39635407224918),c(398.6823159644918,108.89203645203065,308.0042211870933,159.40431574800286),c(315.7349884819466,107.24215711184837,99.16802576099596,195.37246904498812),c(783.9860307221215,672.3258311242802,726.8432946952998,566.4984144275179))
targetgene="Gm10058"
collabel=c("plasmid","Rep1","Rep2","Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(92.31234832734877,127.86564886412691,163.33557184164042,76.02359674135522),c(401.35803620586427,334.92550605700336,812.0111285841551,370.30848735305284),c(458.8860213953715,202.11021917232964,166.83561980967556,123.43616245101761),c(595.3477537053653,184.78648610041566,402.50551632404245,275.48335593372803),c(410.72305705066776,515.5872938069633,691.8428150149483,228.88824825354257),c(325.10000932675007,89.91842403993441,194.83600355395677,188.83280480917264),c(481.6296434470371,607.9805368571712,330.17119165131595,439.79241985859255),c(527.1168875503683,341.5250234177325,124.83504419325374,122.61870442154067),c(591.3341733433067,344.8247820980971,361.6716233636323,442.24479394702337),c(321.0864289646914,426.4938094371201,308.0042211870933,210.08671357557301),c(556.5498102054651,407.52019702502383,287.00393337888244,291.8325165232668),c(433.4666791023334,247.4819010273424,317.33768243518705,223.98350007668097),c(527.1168875503683,311.00225562436026,184.3358596498513,191.28517889760346),c(366.5736730680227,291.2037035421729,161.00220652961696,236.24537051883502),c(354.5329319818468,87.44360502966099,168.00230246568728,65.39664235815502),c(624.780676360462,658.3018567327308,526.1738778612845,476.5780311850547),c(303.6942473957706,166.63781335841054,189.0025902738982,176.57093436701857),c(575.279851895072,335.75044572709453,211.1695607381208,160.2217737774798),c(259.54486341312554,87.44360502966099,225.1697526102614,129.15836865735616),c(181.94897641332514,117.14143315294207,102.66807372903111,72.75376462344747),c(339.8164706542984,256.55623739834493,239.16994448240203,275.48335593372803),c(406.70947668860913,213.6593745536056,376.83849789178464,271.39606578634334),c(464.23746187811633,192.21094313123592,290.50398134691756,222.3485840177271),c(480.29178332635087,395.1461019736567,240.33662713841375,407.094098679515),c(375.9386939128262,95.69300173057239,194.83600355395677,300.8245548475131),c(527.1168875503683,446.29236151930746,325.5044610272691,357.2291588814218))
targetgene="Gm10439"
collabel=c("plasmid","Rep1","Rep2","Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=9
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("Spry2","Vgll4","Nf1","Ndst1","Mga","Bcorl1","Kmt2d","Tead1","Lztr1","Gsk3b")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='Rep1,Rep2,Rep3_vs_plasmid pos.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(461.5617416367439,1851.9895593546123,2279.6979098468955,1619.3843563938137),c(386.6415748783159,2040.9007438054837,2868.8726511328127,2145.0098693474847),c(358.5465123439054,1027.8748289335622,1746.5239360495407,1328.3692979000239),c(596.6856138260515,2716.5263336101284,2572.5352565058365,2377.985407748412))
targetgene="Spry2"
collabel=c("plasmid","Rep1","Rep2","Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(553.8740899640927,1234.9346861264385,1335.8516411334163,795.3866626810604),c(302.3563872750844,612.9301748777181,1068.681312906733,970.3226809891252),c(440.1559797057645,1121.0930116538611,2240.030699542497,1681.511166634061),c(363.8979528266503,1475.8170697930518,899.5123277850339,1221.2822960385452))
targetgene="Vgll4"
collabel=c("plasmid","Rep1","Rep2","Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(155.19177399960085,409.1700763652061,796.8442540560028,348.2371205571755),c(294.3292265509671,1312.479015115006,2185.196614709946,987.4892996081409),c(412.060917171354,2044.2005024858481,1508.5206742231503,1631.6462268359678),c(545.8469292399753,2096.1717017015903,2123.3624339413254,2397.6044004558585))
targetgene="Nf1"
collabel=c("plasmid","Rep1","Rep2","Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(607.3884947915412,1674.627530285017,1796.6912902580445,1411.7500169066716),c(478.9539232056647,1002.3016991607367,999.8470362020416,1519.6544767976275),c(599.361334067424,1333.9274465373755,1456.019954702623,1698.6777852530768),c(468.25104224017497,907.4336371002555,924.01266356128,1159.1554857982978))
targetgene="Ndst1"
collabel=c("plasmid","Rep1","Rep2","Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(765.2559890325145,765.5440138445791,1598.3552387360526,1549.082965858797),c(449.521000550568,938.7813445637188,1012.6805454181705,1419.924597201441),c(339.8164706542984,731.7214873708424,345.3380661794683,753.6963031777367),c(375.9386939128262,839.7885841527818,564.6744055096711,774.13275391466))
targetgene="Mga"
collabel=c("plasmid","Rep1","Rep2","Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(505.71112561938895,1819.8169122210577,1382.5189473738849,1118.282584324451),c(549.8605096020341,1238.2344448068031,1435.019666894412,995.6638799029103),c(389.3172951196883,616.2299335580826,973.0133351137721,993.2115058144794),c(397.34445584380563,666.5512534336422,868.0118960727176,685.0298287016739))
targetgene="Bcorl1"
collabel=c("plasmid","Rep1","Rep2","Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(386.6415748783159,399.27080032411243,876.1786746647996,1395.4008563171328),c(408.0473368092953,1365.275154000839,1258.8505858366427,985.854383549187),c(608.7263549122274,1303.4046787440034,770.0105529677334,1052.068483936819),c(256.86914317175314,479.28994832295314,556.507626917589,550.1492538379791))
targetgene="Kmt2d"
collabel=c("plasmid","Rep1","Rep2","Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(773.2831497566318,1589.6587442656294,1450.1865414225645,1415.0198490245793),c(390.65515524037454,1295.980221713183,1209.8499142841506,1220.4648380090682),c(466.91318211948874,707.7982369381992,798.0109367120145,640.0696370804422),c(480.29178332635087,1544.2870624106165,1276.3508256768187,1173.8697303288827))
targetgene="Tead1"
collabel=c("plasmid","Rep1","Rep2","Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(878.9740992908427,1397.4478011343933,2598.202274938094,1586.6860352147362),c(383.96585463694345,1138.416744725775,974.1800177697838,1021.0050788166953),c(397.34445584380563,1221.7356514049802,1234.3502500603968,1117.465126294974),c(215.3954794304805,103.94239843148381,491.1733981809329,296.73726470012843))
targetgene="Lztr1"
collabel=c("plasmid","Rep1","Rep2","Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(448.18314042988175,681.4001674952827,878.5120399768231,792.9342885926296),c(390.65515524037454,741.6207634119361,1069.8479955627447,931.0846955742321),c(375.9386939128262,1112.0186752828586,634.6753648703741,821.5453196243225),c(568.5905512916411,1046.8484413456583,892.5122318489637,1008.7432083745413))
targetgene="Gsk3b"
collabel=c("plasmid","Rep1","Rep2","Rep3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



dev.off()
Sweave("mouse_output_summary.Rnw");
library(tools);

texi2dvi("mouse_output_summary.tex",pdf=TRUE);

