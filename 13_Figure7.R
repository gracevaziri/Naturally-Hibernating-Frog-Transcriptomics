# A plotting R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# --------------------------------------------------------------------------
# If you don't have the ggplot2 package installed, uncomment the following line:
# install.packages( "ggplot2" );
library( ggplot2 );

# --------------------------------------------------------------------------
# If you don't have the scales package installed, uncomment the following line:
# install.packages( "scales" );
library( scales );

# --------------------------------------------------------------------------
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","plot_X","plot_Y","log_size","value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0000236","mitotic prometaphase",19.735208181154125,-6.918749193006068,2.0166690908116136,1.792391689498254,-2.3788505247882217,1,-0),
                     c("GO:0000278","mitotic cell cycle",2.0818115412710005,1.3108821121214467,-2.8442374015771317,2.5352941200427703,-1.7851142339652717,0.9363876694290494,0.01112148),
                     c("GO:0000280","nuclear division",0.9069880691502312,2.075935628365968,-6.157617456381977,2.1760912590556813,-1.7430324983115053,0.9297529841335284,0.01256073),
                     c("GO:0000819","sister chromatid segregation",0.38349159970781593,1.836420965462401,-5.081550814909317,1.806179973983887,-1.6217896857237633,0.8880342782096778,0.65644738),
                     c("GO:0001817","regulation of cytokine production",0.6695885074263452,-6.928132706041396,-0.49343696920400365,2.0453229787866576,-2.6457351661360327,0.9490618231381225,0.08467301),
                     c("GO:0001906","cell killing",0.15826637448259068,-5.44849699264138,0.48378768839033653,1.4313637641589874,-1.8559298040721304,0.9929565865992958,0.00878426),
                     c("GO:0002250","adaptive immune response",0.31653274896518135,6.076826578157134,2.1013124585426666,1.724275869600789,-1.7430324983115053,0.649369065710227,0.60999175),
                     c("GO:0002252","immune effector process",0.5295836376917458,6.985599012880984,0.97151099068962,1.9444826721501687,-6.631378155520305,0.7180058772507776,0.64525583),
                     c("GO:0002274","myeloid leukocyte activation",0.07304601899196494,6.293306362468327,1.0596594118103384,1.1139433523068367,-1.6217896857237633,0.7175629857042033,0.52780793),
                     c("GO:0002376","immune system process",4.88190893596299,-1.4802081992384155,-0.4459221208998423,2.904715545278681,-1.773399398891063,1,-0),
                     c("GO:0002441","histamine secretion involved in inflammatory response",0.060871682493304116,3.8914379424484244,4.141336405144528,0.3010299956639812,-1.560093910783991,0.6647493434909324,0.58302594),
                     c("GO:0002443","leukocyte mediated immunity",0.17652787923058194,6.760894113554978,0.5609103497078757,1.4771212547196624,-3.4556106965843156,0.6525516965470097,0.57437562),
                     c("GO:0002444","myeloid leukocyte mediated immunity",0.006087168249330412,6.816631316249411,0.22186388635401022,0.3010299956639812,-2.8459133675979653,0.6918704919635426,0.64975923),
                     c("GO:0002819","regulation of adaptive immune response",0.07913318724129534,-2.71849199429194,5.937905830041575,1.146128035678238,-1.307496282507934,0.7044121256601564,0.63603819),
                     c("GO:0002886","regulation of myeloid leukocyte mediated immunity",0.06695885074263452,-1.42940596743662,7.169093517402171,1.0791812460476249,-1.7430324983115053,0.7187224534583754,0.62644919),
                     c("GO:0006428","isoleucyl-tRNA aminoacylation",0.012174336498660824,-4.719793281535931,-3.188469317668655,0.47712125471966244,-1.5245061823833543,0.9646571039744812,0.26938319),
                     c("GO:0006952","defense response",2.672266861456051,3.8847703345253755,5.432985354116361,2.6434526764861874,-3.9174874953910486,0.7605030454510494,0.58409037),
                     c("GO:0006954","inflammatory response",0.9191624056488922,4.229391083355834,5.330193692304753,2.1818435879447726,-5.335234784788577,0.7339863827111892,0.19731568),
                     c("GO:0006955","immune response",2.568785001217434,5.692311844453468,2.4444712566337126,2.6263403673750423,-6.631378155520305,0.6180811331303805,0),
                     c("GO:0009200","deoxyribonucleoside triphosphate metabolic process",0.030435841246652058,-3.6230841762099484,-4.223307470714757,0.7781512503836436,-1.402982359162905,0.9448779529115267,0.54204247),
                     c("GO:0009262","deoxyribonucleotide metabolic process",0.11565619673727782,-3.9398708920105503,-4.358243324664745,1.3010299956639813,-1.8559298040721304,0.940874100304622,0.59958841),
                     c("GO:0009263","deoxyribonucleotide biosynthetic process",0.06695885074263452,-4.100450273653843,-3.7598146559611325,1.0791812460476249,-1.4757114326540062,0.9398790283430708,0.57461369),
                     c("GO:0009607","response to biotic stimulus",2.1852934015096177,3.274394846769108,6.85003545531528,2.5563025007672873,-2.3565787122947492,0.8620526916228055,0.22014536),
                     c("GO:0009611","response to wounding",0.5113221329437546,4.467874135472429,5.748280519562979,1.9294189257142926,-2.246862417089562,0.7922446279219107,0.48716488),
                     c("GO:0015949","nucleobase-containing small molecule interconversion",0.29900292617259927,-4.121915595731434,-4.0951619336743335,2.5440680443502757,-1.909765585185801,0.9392584508067563,-0),
                     c("GO:0030261","chromosome condensation",0.08522035549062576,2.5329140331993862,-5.7779819881402155,1.1760912590556813,-1.4534775588294868,0.9330474672307902,0.67079937),
                     c("GO:0031347","regulation of defense response",1.2235208181154127,-3.619982251220437,5.083386682789907,2.305351369446624,-2.218001560253669,0.7323571004944756,0.55405303),
                     c("GO:0032103","positive regulation of response to external stimulus",0.7730703676649623,-3.463945953170466,4.711418716809515,2.1072099696478683,-1.9261346740906666,0.7127693257320652,0.52459233),
                     c("GO:0034501","protein localization to kinetochore",0.0547845142439737,6.3905307877405315,-3.7110615932395348,1,-1.7430324983115053,0.9637203148030083,0.00808474),
                     c("GO:0038093","Fc receptor signaling pathway",0.024348672997321647,-0.9959549155508418,5.500559718411014,0.6989700043360189,-1.680044471315014,0.6378005816173631,0.57404256),
                     c("GO:0038094","Fc-gamma receptor signaling pathway",0.012174336498660824,-1.2955768987276264,5.502551487448195,0.47712125471966244,-1.6217896857237633,0.6526635678898823,0.54292154),
                     c("GO:0048285","organelle fission",1.0104699293888484,1.6654862947204008,-6.077443337150229,2.2227164711475833,-1.7851142339652717,0.9424315081915543,0.43702484),
                     c("GO:0050729","positive regulation of inflammatory response",0.012174336498660824,-3.968177350973764,4.683484346198842,0.47712125471966244,-2.246862417089562,0.7647511445937802,0.58897718),
                     c("GO:0050776","regulation of immune response",2.3009495982468953,-2.7933308483247,5.642804381792539,2.578639209968072,-3.4556106965843156,0.6284697491413843,-0),
                     c("GO:0051259","protein complex oligomerization",0.5417579741904066,0.7994010509642313,-6.280491566914923,1.954242509439325,-1.356812748664553,0.9587600070353367,0.27281215),
                     c("GO:0060333","type II interferon-mediated signaling pathway",0.048697345994643294,4.786673103465831,3.2455002086368867,1,-1.4300593490027549,0.5606644623392942,0.57246322),
                     c("GO:1903361","protein localization to basolateral plasma membrane",0.024348672997321647,6.051715632097082,-3.9818044873571594,0.6989700043360189,-1.5245061823833543,0.9651715478922231,0.35694904));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$log_size <- as.numeric( as.character(one.data$log_size) );
one.data$value <- as.numeric( as.character(one.data$value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
#head(one.data);


# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = value, size = log_size), alpha = I(0.7) );
p1 <- p1 + scale_colour_gradientn( colours = c(  "#5c004a","#8f0057","#ffa6a6","#f6a772","#ffea58"), limits = c( min(one.data$value), 0) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = log_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) ));
p1 <- p1 + scale_size( range=c(5, 30)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.15, ];
# p1 <- p1 + geom_text( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 5 , fill = "white");
p1 <- p1 + geom_label(data = ex, aes(plot_X, plot_Y, label = description),fill = "white", colour = I(alpha("black", 0.85)), size = 5 );

p1 <- p1 + labs (y = "semantic space x", x = "semantic space y",size = "Generality of GO term",
                 color = "-log10 B-H \nadjusted p value ");
p1 <- p1 + theme(legend.key = element_blank(), panel.grid = element_blank(), axis.title =  element_text(size = 16),
                 axis.text =  element_text(size = 14)) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/5,max(one.data$plot_X)+one.x_range/10);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);


p1
# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).

