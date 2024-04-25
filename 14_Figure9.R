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
revigo.data <- rbind(c("GO:0001934","positive regulation of protein phosphorylation",0.3713172632091551,-3.6940504969339685,-4.303003573639683,1.792391689498254,-1.6840997137818319,0.7122381749206305,0.68828314),
                     c("GO:0001991","regulation of systemic arterial blood pressure by circulatory renin-angiotensin",0.006087168249330412,3.5728883508952154,5.27057119203892,0.3010299956639812,-1.3760896480546199,0.9347392029693122,0.18823519),
                     c("GO:0002011","morphogenesis of an epithelial sheet",0.06695885074263452,5.763075753022557,3.2528902965355733,1.0791812460476249,-2.05784738820196,0.9685037854835181,-0),
                     c("GO:0002376","immune system process",4.88190893596299,7.2544268801118355,1.8327389443034479,2.904715545278681,-1.9378392788576968,1,-0),
                     c("GO:0002523","leukocyte migration involved in inflammatory response",0.018261504747991236,4.427201058427495,-6.179291045319503,0.6020599913279624,-2.444483438508703,0.8894769490093627,0.55871077),
                     c("GO:0002682","regulation of immune system process",2.81227173119065,-0.5887553850176576,0.9758245725498381,2.6655809910179533,-1.3621227405498537,0.9082293730919837,0.0906811),
                     c("GO:0002920","regulation of humoral immune response",0.1643535427319211,-0.4952621534395268,-8.02947728557454,1.4471580313422192,-2.442662458714289,0.637880464344814,0.6816911),
                     c("GO:0006508","proteolysis",5.575846116386657,-1.2789265787440027,5.034286931797268,2.962369335670021,-4.202379470532863,0.9615137362409483,0.33044758),
                     c("GO:0006534","cysteine metabolic process",0.018261504747991236,-2.416805207641587,5.458302792954209,0.6020599913279624,-1.409588037526385,0.96866334941829,0.47179469),
                     c("GO:0006952","defense response",2.672266861456051,5.333280560630943,-4.275032188736024,2.6434526764861874,-2.797600785113235,0.9113857602641301,0.54559397),
                     c("GO:0006954","inflammatory response",0.9191624056488922,5.175656163503305,-4.714176976013962,2.1818435879447726,-1.4889072262528242,0.9123711176541058,0.58409037),
                     c("GO:0006955","immune response",2.568785001217434,4.257807878436368,-5.1819663337487,2.6263403673750423,-2.100094725438057,0.8731797422908455,0.6035129),
                     c("GO:0006957","complement activation, alternative pathway",0.018261504747991236,1.2775794155221083,-5.450390066469358,0.6020599913279624,-4.284721443615654,0.53563023782914,-0),
                     c("GO:0007039","protein catabolic process in the vacuole",0.006087168249330412,-0.17484665675376418,5.864469008038238,0.3010299956639812,-1.6960232387441632,0.9700419599499517,0.31880935),
                     c("GO:0007586","digestion",0.1643535427319211,3.336669686846259,4.5521088319836505,1.4471580313422192,-16.615473332381555,0.986359577037019,0),
                     c("GO:0009605","response to external stimulus",3.627952276600925,6.173571709067233,-3.2473044333677263,2.775974331129369,-2.2146478812703045,0.935148555756214,0.13818646),
                     c("GO:0009611","response to wounding",0.5113221329437546,5.537999543606507,-4.842078535172162,1.9294189257142926,-3.2376470916805107,0.9228238800696803,0.34966636),
                     c("GO:0009612","response to mechanical stimulus",0.1095690284879474,3.6071246605680325,-2.007885446306482,1.2787536009528289,-2.2839930650987568,0.9156765250014612,0.41442906),
                     c("GO:0009617","response to bacterium",0.6939371804236669,4.245230367070538,-2.644052046680374,2.060697840353612,-1.5192659495652874,0.9017297326941557,0.62014566),
                     c("GO:0009725","response to hormone",1.065254443632822,6.5024098783010995,-2.0704056661427974,2.24551266781415,-1.3376160490506184,0.9174514227114441,0.67278961),
                     c("GO:0009991","response to extracellular stimulus",0.4808862916971025,4.240373912784125,-2.2615804871685103,1.9030899869919435,-1.6865447324991856,0.9060332319419515,0.59740402),
                     c("GO:0010466","negative regulation of peptidase activity",0.012174336498660824,-5.488921977776669,-0.6192366011114309,0.47712125471966244,-1.9454574953032397,0.8065438101335516,0.42665476),
                     c("GO:0010574","regulation of vascular endothelial growth factor production",0.018261504747991236,-5.91956325806181,-2.903584166509165,0.6020599913279624,-1.6582314247110288,0.8240337543125656,0.43752337),
                     c("GO:0010575","positive regulation of vascular endothelial growth factor production",0.018261504747991236,-4.850837298003668,-3.8663663028065227,0.6020599913279624,-1.7432599627349825,0.7242245631894424,0.67592971),
                     c("GO:0010758","regulation of macrophage chemotaxis",0.012174336498660824,-1.0911685062790375,-7.781101507826363,0.47712125471966244,-2.9242851925069995,0.6328150858150219,0.39391674),
                     c("GO:0010866","regulation of triglyceride biosynthetic process",0.012174336498660824,-3.733610045815214,-0.7467877401710794,0.47712125471966244,-1.680074948720387,0.900495557534933,0.04698262),
                     c("GO:0014061","regulation of norepinephrine secretion",0.006087168249330412,-6.604869989947729,2.1815495109576997,0.3010299956639812,-2.444483438508703,0.8577924879684902,0.6878551),
                     c("GO:0014070","response to organic cyclic compound",0.5356708059410762,6.720975854761098,-1.81266134133374,1.9493900066449128,-1.8066288146690812,0.9217751686894737,0.62718907),
                     c("GO:0019835","cytolysis",0.0547845142439737,2.128011168663172,1.895167923267097,1,-1.3331782506057928,0.9989133653563766,-0),
                     c("GO:0030574","collagen catabolic process",0.18870221572924276,0.9193894364339499,6.466099471088072,1.505149978319906,-2.2210420734326797,0.9788263104607441,-0),
                     c("GO:0031214","biomineral tissue development",0.09739469198928659,5.745649484927123,3.4878177170094036,1.2304489213782739,-1.680074948720387,0.9732384196975973,0.51393908),
                     c("GO:0031667","response to nutrient levels",0.46871195519844167,3.9460001473952513,-2.5372004034506195,1.8920946026904804,-1.6556799509726938,0.8998413099127641,0.52485425),
                     c("GO:0032101","regulation of response to external stimulus",1.5826637448259069,-0.16044591633076172,-8.043796896155246,2.416640507338281,-1.5603090425904729,0.756128339039395,0.35514145),
                     c("GO:0032637","interleukin-8 production",0.24174753904483634,1.2494973571301229,4.6287190344741145,3.1405080430381798,-2.2146478812703045,0.9570138999134599,0.243963),
                     c("GO:0032642","regulation of chemokine production",0.018261504747991236,-5.831131816827102,-2.712286034087643,0.6020599913279624,-1.4155602155814344,0.8240337543125656,0.57767955),
                     c("GO:0032722","positive regulation of chemokine production",0.018261504747991236,-4.769533436540296,-3.7450666399197488,0.6020599913279624,-2.1667190490768347,0.7242245631894424,0.58148276),
                     c("GO:0032963","collagen metabolic process",0.23739956172388604,-4.385454954194872,3.302387699773671,1.6020599913279623,-1.8507186429296238,0.9836612596999027,0.09567583),
                     c("GO:0033494","ferulate metabolic process",0.34702264427065216,-2.9113845116637815,5.1452164961519395,1.2304489213782739,-1.6138512687089537,0.9696996274519901,0.16944045),
                     c("GO:0033602","negative regulation of dopamine secretion",0.012174336498660824,-6.6045877699559306,1.4570493656611254,0.47712125471966244,-2.9905533536788664,0.8061790412346083,0.04217913),
                     c("GO:0033993","response to lipid",0.5234964694424153,6.5825701720011365,-1.597150772630833,1.9395192526186185,-2.1352265822757164,0.9219143363597697,0.19230277),
                     c("GO:0040012","regulation of locomotion",0.9800340881421963,-5.6580759496952036,4.038502539637359,2.2095150145426308,-1.9898999751051603,0.9164160008457223,0.05616585),
                     c("GO:0040017","positive regulation of locomotion",0.5174093011930849,-3.3754599988630294,-6.085988245426478,1.9344984512435677,-2.6372542763206925,0.7082725740226002,0.64156167),
                     c("GO:0045745","positive regulation of G protein-coupled receptor signaling pathway",0.024348672997321647,-2.2611923416625506,-6.52758923670197,0.6989700043360189,-1.680074948720387,0.7341057628511282,0.37653755),
                     c("GO:0045766","positive regulation of angiogenesis",0.09130752373995618,-4.975471444985081,-4.5007181705415285,1.2041199826559248,-2.1667190490768347,0.7366109398032381,0.40875019),
                     c("GO:0045917","positive regulation of complement activation",0.04969119841393834,-1.9341220461802022,-6.7262810021665,1.0413926851582251,-2.9905533536788664,0.5526535558556124,0.49601603),
                     c("GO:0048521","negative regulation of behavior",0.006087168249330412,-6.58304397619869,-1.8042254317789939,0.3010299956639812,-2.2839930650987568,0.8270247929649193,0.23108699),
                     c("GO:0048709","oligodendrocyte differentiation",0.0547845142439737,4.584401505429577,4.053701066496406,1,-1.6865447324991856,0.9695732174348984,0.26599572),
                     c("GO:0050776","regulation of immune response",2.3009495982468953,-0.42470108731100686,-7.813240050198303,2.578639209968072,-2.3135228224615885,0.5714604279092712,0.6100295),
                     c("GO:0050795","regulation of behavior",0.0547845142439737,-6.749805278591267,-3.0656415404591755,1,-1.6582314247110288,0.8542611426488598,0.47065417),
                     c("GO:0050900","leukocyte migration",0.28609690771852936,3.368725816248366,-7.640015616208324,1.6812412373755872,-2.5133400362534646,0.9263032826410335,0.39713414),
                     c("GO:0072376","protein activation cascade",0.030435841246652058,-0.4844501408596177,4.713408651871965,0.7781512503836436,-1.9378392788576968,0.9689600091031542,0.32119648),
                     c("GO:1901615","organic hydroxy compound metabolic process",1.467007548088629,-2.0731249595430583,3.8667366292831415,2.383815365980431,-1.5941392296388364,0.9752628009577852,0.11479035),
                     c("GO:2000259","positive regulation of protein activation cascade",0.05847728127063996,-3.705268499082215,-4.056081481716543,0.47712125471966244,-2.9905533536788664,0.7082884359816564,0.24220904));

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
ex <- one.data [ c(one.data$dispensability < 0.05), ];
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

