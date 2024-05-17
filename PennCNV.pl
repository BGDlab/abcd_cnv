use warnings;
#use strict;
@ARGV>=1 or die "Usage: $0 <input-file>\n the input file contains one int.csv path per line
OPTIONS:
--keepSignalFile
--clean_cnv
--chrx
--chry
--build
--GCAdjustOut
--noGC
--intersectionPFB
--noPennCNV
--loh
--minLenLoh
--conf
--qcFilterPFB
";
use Getopt::Long;
GetOptions('keepSignalFile'=>\ my $keepSignalFile,'clean_cnv'=>\ my $clean_cnv,'chrx'=>\ my $chrx,'chry'=>\ my $chry,'build=s'=>\ (my $build = 'hg19'),'GCAdjustOut'=>\ my $GCAdjustOut,'noGC'=>\ my $noGC,'intersectionPFB'=>\ my $intersectionPFB,'noPennCNV'=>\ my $noPennCNV,'loh'=>\ my $loh,'minLenLoh=i'=>\ (my $minLenLoh = 0),'conf'=>\ my $conf, 'qcFilterPFB'=>\ my $qcFilterPFB);
open (FAM, $ARGV[0]) or die;
open (LOG, ">$ARGV[0].Intersect_RunPennCNV.log");
unless(-d 'temp')
{
       mkdir 'temp' or die;
}
while (<FAM>) {
	s/[\r\n]+$//;
	my @ind = split (/\s+/, $_);
	map {s/\.ped/.int.csv/} @ind;
	#my @triofile;
	my $i=0;
	$ind[$i] =~ m/([^\\\/]+)$/;
	#push @triofile, $1;
	my $command;
	my @out;
	#mv or cp 2 column .int.csv in cwd into temp, then make 5 column .int.csv in cwd
	$cmd="cp /mnt/isilon/cag_ngs/hiseq/glessner/ABCD_Array_ForCNVCalling_Aaron/genomics_sample03/CNV/ABCD_Array_ForCNVCalling_LRR_BAF_int.csv_x11088_5col/$ind[$i] $ind[$i]";
	@out=`$cmd`;
	#Zip not working on HPC9 Respublica, untarring specific files is much slower
	#$cmd="tar -zxvf ABCD_Array_ForCNVCalling_LRR_BAF_int.csv_x11088.tar.gz $ind[$i] > temp/$ind[$i]";
	#@out=`$cmd`;
	#$cmd="unzip -p /mnt/isilon/cag_inf_op/glessner/CAG_Internal_IntCsv_BafLrr_SignalFiles_Zip/AllEpicPatIDsSIDs_Paths_98CR_IntCsv_2colIntCsvs_ForRespublica.zip $ind[$i]>temp/$ind[$i]";
	#@out=`$cmd`;
=pod	
my $MatchFirstZip=`fgrep $ind[$i] ../AllEpicPatIDsSIDs_Paths_98CR_IntCsv_2colIntCsvs_ForRespublica.zip_ChipIDs.int.csv.txt | wc -l`;
        if($MatchFirstZip==1)
        {
		print "Unzip 1\n";
                $cmd="unzip -p /mnt/isilon/cag_inf_op/glessner/CAG_Internal_IntCsv_BafLrr_SignalFiles_Zip/AllEpicPatIDsSIDs_Paths_98CR_IntCsv_2colIntCsvs_ForRespublica.zip $ind[$i] > temp/$ind[$i]";
        }
        else
        {
		print "Unzip 2\n";
                $cmd="unzip -p /mnt/isilon/cag_inf_op/glessner/Munir_NDD_CNV_DeepCNV_IntCsv_2colIntCsvs_ForRespublica/Munir_NDD_CNV_DeepCNV_IntCsv_2colIntCsvs_ForRespublica.zip $ind[$i] > temp/$ind[$i]";
        }
        @out=`$cmd`;
        if(@out ne ""){print "OUT UNZIP:".@out;}
	#$ind[$i]="temp/".$ind[$i];
=cut
	#open(OUT,">$ind[$i]"); print OUT "Name\tChr\tPosition\tsample.Log R Ratio\tsample.B Allele Freq\n"; close (OUT);
		$command = "wc -l $ind[$i]";
		print LOG "running $command\n";
		my $SnpCount= `$command`;
		print LOG "SnpCount:$SnpCount";
		my @CountPath = split (/\s+/, $SnpCount);
        #$build = $build eq "hg19" ? "hg19" : "hg18";
        my %platver = (
        "BDCHP-1X10-HUMANHAP550", 555352,
        "HumanHap550v3", 561466,
        "Human610-Quadv1", 620901,
        "HumanOmniExpress-12v1", 733202,
        "Human660W-Quad_v1", 657366,
        "HumanOmni1-Quadv1-0", 1134514,
        "HumanOmni25-8v1", 2379855,
        "HumanOmni25-4v1", 2443179,
        "Human1M-Duov3-0", 1199187,
        "HumanOmniExpressExome-8v1", 951117,
        "IBD_Consort_ExomePlus", 270495,
        "HumanExome-12v1", 247870,
        "Cardio-Metabo_Chip", 196725,
        "CVDSNP55", 49094,
        "LungCa_GxE", 11376,
		"HumanOmniExpress-12v1-1", 719665,
		"HumanOmniExpress-24v1-0", 716503,
        "CAG_AUT_86K", 8582,
		"iSelect_CVD3",	53831,
		"HumanOmni25-8v1-1",	2391739,
		"HumanExome-12v1-1",	242901,
		"CVD_Upenn_Phase1_49234",	45707,
		"HumanOmniExpressExome-8v1-2",	964193,
		"hakonhakonarson35277714",	22459,
		"CytoSNP-850K",	851622,
		"HumanCoreExome-12v1-0",	538448,
		"NeuroX",	267607,
		"CAG_PSP_5960",	5283,
		"HumanCore-12v1-0",	298930,
		#"CytoConsortiumArray",	851622,
		"CanineHD",	173662,
		"HumanExon510S-2v1.0",	511354,
		"HumanCytoSNP-12v1-0", 299671,
		"HumanCNV-12_v1-0",	52167,
		"InfiniumOmni25Exome-8v1-2", 2608742,
		"InfiniumOmni2-5-8v1-3", 2372784,
		"HumanOmniExpress-24v1-1", 713014,
		"GSAMD-24v1-0", 700078,
		"HumanOmni25-8v1-2", 2338671,
		"HumanCore-24v1-0", 306670,
		"HumanCytoSNP_12v2_1_L", 294602,
		"HumanCytoSNP-12v2",	301232,
		"Human1M-Single",1072820,
		"CytoSNP-850Kv1-1_iScan_B1",843888,
		"Immuno_BeadChip_11419691_B",196524,
		"Human1M-Duov3_H",1192666,
		"HumanOmni2.5-4v1_H",2443177,
		"HumanCoreExome-24v1-0_A",547644,
		"InfiniumOmniExpressExome-8v1-4_A1",960919,
		"InfiniumOmniExpressExome-8v1-3_A",958497,
		"PsychArray-B",571054,
		"CytoSNP-850K_B",851275,
		"Exomev1-1-Neurov1-0-NeuroX",449049,
		"MEGA_Consortium_v2_15070954_A2",2036060,
		"GSA-mGluR_enrichd_20011739X343186_B2",683369,
		"GSAMD-24v2-0",759993,
		"CytoSNP-850Kv1-2",846500,
		"Multi-EthnicGlobal",1779819,
		"NateraSMARTARRAY",655089,
		"Omni2-5-8v1-4",2382209,
		"DragonfishBovineBeta96",7917,
		"UNIVERSITY-HOMOSAPIEN-81",312960,
		"OmniExpressExome-8v1-6",962215,
		"GDA-8v1-0",1914935,
		"GSA-24v2-0",665608,
		"ABCD_Affymetrix_SmokeScreen_BATCH_3_Saliva.sample",611122,
"ABCD_Affymetrix_SmokeScreen_BATCH_4_Saliva.sample",611293,
"ABCD_Affymetrix_SmokeScreen_BATCH_2_Saliva.sample",611359,
"ABCD_Affymetrix_SmokeScreen_BATCH_1_Saliva.sample",614706,
"ABCD_Affymetrix_SmokeScreen_BATCH_2_WB.sample",664231,
		"GSA-24v3",654027
        ) ;
	my $fg=0;
    foreach my $key (keys %platver) {
        if($CountPath[0]-1 == $platver{$key}) {
            #$command = qq(awk \'\{print \$2"\\t\"\$1\"\\t\"\$4\}\' \/mnt\/isilon\/cag_ngs\/hiseq\/glessner\/mnt\/public\/SNPMAP_FILES\/CLEAN_MAPs\/$key\_$build.map \> NameChrPos.txt);
            #`$command`;
            #$command = "paste NameChrPos.txt temp/$ind[$i] | perl -pe 's/,/\t/' | sed 's/\r//' | grep -v REMOVE >> $ind[$i]";
            #`$command`;
	$fg=$key;	
        }
    }
		if ( $fg eq "0")
		{
			print LOG "ERROR: SnpCount $CountPath[0] definition files not found!\n";
		}
	
		my $PFB_FILE="/mnt/isilon/cag_ngs/hiseq/glessner/mnt/public/SNPMAP_FILES/PFB\_".$build."/".$fg."\_".$build.".pfb";
		if($qcFilterPFB)
		{
		$cmd="fgrep -vf /mnt/isilon/cag_ngs/hiseq/glessner/mnt/public/SNPMAP_FILES/QC/$fg.lmiss.1 $PFB_FILE > temp.pfb";
		print LOG "running $cmd\n";
		$out=`$cmd`;
		print "$out\n";
		$PFB_FILE="temp.pfb";
		}
			
		if(!($noPennCNV))
		{
		$command = "perl /cm/shared/easybuild/software/PennCNV/1.0.5-GCCcore-11.3.0/detect_cnv.pl -test -hmm /cm/shared/easybuild/software/PennCNV/1.0.5-GCCcore-11.3.0/lib/hhall.hmm -pfb $PFB_FILE --gcmodel /mnt/isilon/cag_ngs/hiseq/glessner/mnt/public/SNPMAP_FILES/GC\_$build/$fg\_$build.gcmodel --conf $ind[$i] >myCall.rawcnv 2>>$ARGV[0].Intersect.log";             #BUILD SPECIFICATION AND SNP SET TO INCLUDE IN PFB
		print LOG "running $command\n";
		my @output=`$command`;
		if($chrx)
		{
		$command = "perl /cm/shared/easybuild/software/PennCNV/1.0.5-GCCcore-11.3.0/detect_cnv.pl -test -hmm /cm/shared/easybuild/software/PennCNV/1.0.5-GCCcore-11.3.0/lib/hhall.hmm -pfb $PFB_FILE --gcmodel /mnt/isilon/cag_ngs/hiseq/glessner/mnt/public/SNPMAP_FILES/GC\_$build/$fg\_$build.gcmodel --chrx --conf $ind[$i] >>myCall.rawcnv 2>>$ARGV[0].Intersect.log";             #BUILD SPECIFICATION AND SNP SET TO INCLUDE IN PFB
		print LOG "running $command\n";
		my @output=`$command`;
		print "PennCNV X output:"."@output\n";
		}
		if($chry)
		{
		$command = "awk '{if(\$2==\"Y\")print \$1\"\\tX\\t\"\$3\"\\t\"\$4}' $PFB_FILE > myYasX.pfb";
		print LOG "running $command\n";
                @output=`$command`;
		print "PennCNV Y output:"."@output\n";
		$command = "perl /cm/shared/easybuild/software/PennCNV/1.0.5-GCCcore-11.3.0/detect_cnv.pl -test -hmm /cm/shared/easybuild/software/PennCNV/1.0.5-GCCcore-11.3.0/lib/hhall.hmm -pfb myYasX.pfb --gcmodel /mnt/isilon/cag_ngs/hiseq/glessner/mnt/public/SNPMAP_FILES/GC\_$build/$fg\_$build.gcmodel --chrx --conf $ind[$i] >myYasXCall.rawcnv 2>myYasXCall.log";
		print LOG "running $command\n";
                @output=`$command`;
		$command = "sed 's/X/Y/' myYasXCall.rawcnv >>myCall.rawcnv";
		print LOG "running $command\n";
                @output=`$command`;
		$command = "sed 's/X/Y/g' myYasXCall.log >>$ARGV[0].Intersect.log";
		print LOG "running $command\n";
                @output=`$command`;
		}
		
		if($clean_cnv)
		{
		#$command = "cat  myCall.rawcnv >>$ARGV[0].Intersect.log";
		#print LOG "running $command\n";
                #@output=`$command`;
		$command = "perl /cm/shared/easybuild/software/PennCNV/1.0.5-GCCcore-11.3.0/clean_cnv.pl combineseg myCall.rawcnv --signalfile $ind[$i] >> $ARGV[0].Intersect_Clean_CNV.rawcnv 2>>$ARGV[0].Intersect.log";
		print LOG "running $command\n";
		@output=`$command`;
		}
		
		if($loh)
		{
		$command = "perl /cm/shared/easybuild/software/PennCNV/1.0.5-GCCcore-11.3.0/detect_cnv.pl -test -hmm /cm/shared/easybuild/software/PennCNV/1.0.5-GCCcore-11.3.0/lib/hh550.hmm -pfb $PFB_FILE --gcmodel /mnt/isilon/cag_ngs/hiseq/glessner/mnt/public/SNPMAP_FILES/GC\_$build/$fg\_$build.gcmodel --loh --conf $ind[$i] >>myCall.rawcnv 2>>$ARGV[0].Intersect.log";
                print LOG "running $command\n";
                @output=`$command`;
		}
		$command = "cat myCall.rawcnv >> $ARGV[0]_Intersect.rawcnv";
		print LOG "running $command\n";
		@output=`$command`;
		$command = "grep cn=2 myCall.rawcnv | awk '{if(NF<8)print}' | sed 's/length=//' | awk '{gsub(/,/,\"\",\$3); if(\$3 > $minLenLoh)print \$1,\$2,\"length=\"\$3,\$4,\$5,\$6,\$7}' >> $ARGV[0].Intersect_Clean_CNV.rawcnv";
		print LOG "running $command\n";
		@output=`$command`;
		
	
		}
		
		if($GCAdjustOut)
		{
		$command = "perl /cm/shared/easybuild/software/PennCNV/1.0.5-GCCcore-11.3.0/genomic_wave.pl -adjust -gcmodel /mnt/isilon/cag_ngs/hiseq/glessner/mnt/public/SNPMAP_FILES/GC\_$build/$fg\_$build.gcmodel $ind[$i]";
		print LOG "running $command\n";
		@output=`$command`;
		print @output;
		}
		
		if(!($keepSignalFile))
		{
			$command = "rm $ind[$i]";
			print LOG "running $command\n";
			`$command`;
		}
		
}
		if($conf)
		{
		$command = "awk '{print \$1,\$2,\$3,\$4,\"./\"\$5,\$6,\$7,\$8}' $ARGV[0].Intersect_Clean_CNV.rawcnv | sed 's/chrX/chr23/' | sed 's/chrY/chr24/' | sed 's/state4,cn=2/state5,cn=3/' > $ARGV[0].Intersect_Clean_CNV.rawcnv_wPath";
		print LOG "running $command\n";
                @output=`$command`;
		
		$command = "awk '{print \$1,\$2,\$3,\$4,\"./\"\$5,\$6,\$7,\$8,\$9,\$10,\$11,\$12,\$13,\$14}' $ARGV[0].Intersect.log > $ARGV[0].Intersect.log_wPath";
		print LOG "running $command\n";
                @output=`$command`;
		$command = "perl /cm/shared/easybuild/software/PennCNV/1.0.5-GCCcore-11.3.0/filter_cnv.pl $ARGV[0].Intersect_Clean_CNV.rawcnv_wPath -qclogfile $ARGV[0].Intersect.log_wPath -qcsumout QCsum.qcsum -out goodCNV.good.cnv -chroms 1-22";
		print LOG "running $command\n";
                @output=`$command`;
		$command = "Rscript /work/libspace1a1/penncnv2011Jun16/R_script_convert_raw_cnv.R $ARGV[0].Intersect_Clean_CNV.rawcnv_wPath QCsum.qcsum ./R";
		print LOG "running $command\n";
                @output=`$command`;
		$command = "mkdir R";
		print LOG "running $command\n";
                @output=`$command`;

		$command = "Rscript /work/libspace1a1/penncnv2011Jun16/R_script_calculate_quality_score.R ./R";
		print LOG "running $command\n";
                @output=`$command`;

		$command = "rm $ARGV[0].Intersect_Clean_CNV.rawcnv_wPath $ARGV[0].Intersect.log_wPath QCsum.qcsum goodCNV.good.cnv myCall.rawcnv myYasXCall.log myYasXCall.rawcnv myYasX.pfb NameChrPos.txt";
		print LOG "running $command\n";
                @output=`$command`;

		}	

