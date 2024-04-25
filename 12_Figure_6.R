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
revigo.data <- rbind(c("GO:0001773","myeloid dendritic cell activation",0.006087168249330412,5.374311906862038,0.7617796052152965,0.3010299956639812,-2.6210990450269445,0.7160448266035724,0.66562363),
                     c("GO:0002250","adaptive immune response",0.31653274896518135,4.4618430311064055,2.7995666875030674,1.724275869600789,-13.070407296448792,0.692174278183734,0.60999175),
                     c("GO:0002252","immune effector process",0.5295836376917458,4.918796788107134,1.9118811408861944,1.9444826721501687,-21.90001004864137,0.7054513270720102,0.64525583),
                     c("GO:0002262","myeloid cell homeostasis",0.14000486973459947,5.465819868183063,2.3143482134435023,1.380211241711606,-2.981513373300524,0.7120416887375761,0.56136464),
                     c("GO:0002318","myeloid progenitor cell differentiation",0.012174336498660824,1.2423202245341711,-5.871677142082978,0.47712125471966244,-1.9231348588820418,0.8563895834030827,0.68785068),
                     c("GO:0002335","mature B cell differentiation",0.006087168249330412,3.7420056720896056,-1.2474887602925522,0.3010299956639812,-5.3915688605487695,0.632604886938222,0.66562363),
                     c("GO:0002437","inflammatory response to antigenic stimulus",0.012174336498660824,4.726532687749406,4.085782296883788,0.47712125471966244,-7.84117148706337,0.7428992091716821,0.57237271),
                     c("GO:0002438","acute inflammatory response to antigenic stimulus",0.07304601899196494,4.7188735103760155,3.6364612167570916,1.255272505103306,-9.645504944978146,0.6976821137707256,0.52780793),
                     c("GO:0002443","leukocyte mediated immunity",0.17652787923058194,5.5416195419440655,1.749565432172194,1.4771212547196624,-12.592471226382264,0.6832011165145655,0.57437562),
                     c("GO:0002444","myeloid leukocyte mediated immunity",0.006087168249330412,6.0712046850662755,1.3476715838165805,0.3010299956639812,-10.499311131461612,0.732218658374162,0.64975923),
                     c("GO:0002448","mast cell mediated immunity",0.006087168249330412,6.170556864379345,1.6544876527086525,0.3010299956639812,-1.5518741833119813,0.7241796445675948,0.65194433),
                     c("GO:0002460","adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains",0.18261504747991236,4.725313731787601,2.9689250988766047,1.4913616938342726,-11.156253939905444,0.6903887883683634,0.57632924),
                     c("GO:0002520","immune system development",0.14000486973459947,4.063891679590744,-0.12286528258709797,1.380211241711606,-14.861719987165944,0.6913813909800086,0.56136464),
                     c("GO:0002521","leukocyte differentiation",0.23739956172388604,1.3588364154989203,-5.6709445316348095,1.6020599913279623,-12.330350463777501,0.8308544379503954,0.55325109),
                     c("GO:0002579","positive regulation of antigen processing and presentation",0.19997508619561286,-4.838807478265908,3.160988106801779,1.8129133566428555,-1.4126177530671056,0.4536563321460942,0.68472405),
                     c("GO:0002682","regulation of immune system process",2.81227173119065,-4.309012763452796,-4.851891458080947,2.6655809910179533,-38.294054661351254,0.9166876973328667,-0),
                     c("GO:0002684","positive regulation of immune system process",1.9783296810323838,-4.678765719699628,2.961169365035176,2.513217600067939,-36.26196333535009,0.38832950696893925,0.09918063),
                     c("GO:0002685","regulation of leukocyte migration",0.14000486973459947,-4.933644748191695,2.2746674038673,1.380211241711606,-5.457121326074864,0.47579275642070085,0.6614621),
                     c("GO:0002688","regulation of leukocyte chemotaxis",0.060871682493304116,-4.970243372310637,2.008332724704772,1.0413926851582251,-4.623569638813699,0.480842564342099,0.67298281),
                     c("GO:0002699","positive regulation of immune effector process",0.20696372047723396,-4.677144389143841,3.349583106643187,1.5440680443502757,-17.08682218935012,0.39617923152995554,0.68705208),
                     c("GO:0002819","regulation of adaptive immune response",0.07913318724129534,-4.571287081232337,3.9812838813114277,1.146128035678238,-12.89138760912091,0.4927048448411847,0.62735422),
                     c("GO:0002861","regulation of inflammatory response to antigenic stimulus",0.006087168249330412,-5.0279948956833405,1.5730153410951333,0.3010299956639812,-7.4790917171685765,0.5366802699959387,0.50929318),
                     c("GO:0006955","immune response",2.568785001217434,4.422632948005342,2.1995572294656514,2.6263403673750423,-49.31454451488888,0.6630185915123157,0),
                     c("GO:0019882","antigen processing and presentation",0.08522035549062576,5.7015844155012045,2.673925435228065,1.1760912590556813,-13.151452142047942,0.7399869704453812,0.535391),
                     c("GO:0030097","hemopoiesis",0.547845142439737,1.8504527492577256,-5.849083036846654,1.9590413923210936,-15.865203864993562,0.8846936702199997,-0),
                     c("GO:0030219","megakaryocyte differentiation",0.012174336498660824,0.935657368691419,-5.757197375642068,0.47712125471966244,-2.254206795626099,0.8506887237743037,0.68785068),
                     c("GO:0030225","macrophage differentiation",0.006087168249330412,0.5304874126199098,-5.654290919243429,0.3010299956639812,-3.162085962987346,0.8471292045959163,0.65757868),
                     c("GO:0030595","leukocyte chemotaxis",0.15217920623326028,4.078230155508558,2.47858417226141,1.414973347970818,-8.41823624864301,0.7049760405258867,0.65467397),
                     c("GO:0034340","response to type I interferon",0.06695885074263452,4.444457558787007,3.410173859890073,1.0791812460476249,-2.067602801694835,0.7151734077313723,0.61642722),
                     c("GO:0034341","response to type II interferon",0.0547845142439737,4.970101506661381,3.436873127000759,1,-9.965685611998142,0.7186949456501772,0.6034586),
                     c("GO:0035855","megakaryocyte development",0.012174336498660824,0.7429161410788532,-5.709317126827291,0.47712125471966244,-2.2648539925601106,0.8497318122527764,0.68785068),
                     c("GO:0038093","Fc receptor signaling pathway",0.024348672997321647,-3.6859813611847896,3.8709201433742826,0.6989700043360189,-9.574811564515493,0.4946635609596797,0.56695953),
                     c("GO:0038094","Fc-gamma receptor signaling pathway",0.012174336498660824,-3.762322939324535,4.2734416508178255,0.47712125471966244,-10.022056269888916,0.5127288613826825,0.53658145),
                     c("GO:0043374","CD8-positive, alpha-beta T cell differentiation",0.012174336498660824,3.767078574712164,-1.0773192903717947,0.47712125471966244,-1.5111002115787675,0.6193860018896508,0.69806793),
                     c("GO:0045321","leukocyte activation",0.547845142439737,4.896866868264447,1.4914992699802538,1.9590413923210936,-20.618942377017316,0.6383423893609429,0.6477224),
                     c("GO:0048002","antigen processing and presentation of peptide antigen",0.048697345994643294,5.6228608119489705,3.2064938202122035,0.9542425094393249,-11.595904209290216,0.7203753479839076,0.50885081),
                     c("GO:0048534","hematopoietic or lymphoid organ development",0.03652300949598247,4.28377018615074,-0.39715679128554493,0.8450980400142568,-15.939870785956476,0.710735203332712,0.49620584),
                     c("GO:0050853","B cell receptor signaling pathway",0.09739469198928659,0.703788938357521,2.8668183832999636,1.2304489213782739,-7.732557305180571,0.3088368631784817,0.63935224),
                     c("GO:0050854","regulation of antigen receptor-mediated signaling pathway",0.07304601899196494,-4.676227142780233,3.7079475073839374,1.1139433523068367,-2.884446066771061,0.4949592703403742,0.6228485),
                     c("GO:0050857","positive regulation of antigen receptor-mediated signaling pathway",0.030435841246652058,-5.035596174449777,2.774217621478941,0.7781512503836436,-1.5688679696336858,0.5023213379037502,0.57748458),
                     c("GO:0050900","leukocyte migration",0.28609690771852936,4.943514291742487,2.259025520927106,1.6812412373755872,-14.29953594613169,0.7061670586103083,0.6035129),
                     c("GO:0060333","type II interferon-mediated signaling pathway",0.048697345994643294,3.6139796479803503,3.396925641182087,1,-10.765921222667215,0.6596625001479184,0.50885081),
                     c("GO:1902105","regulation of leukocyte differentiation",0.15217920623326028,-4.7673147252415715,2.553159389821536,1.414973347970818,-14.782610078036761,0.4224796903677156,0.66675985));

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

