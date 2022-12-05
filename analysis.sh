BASEDIR=/media/ggj/Files/mount/P53/NvP53

# CodeTest --------------------
MOTIFDB=$BASEDIR/Resources/0_TFBS/JASPAR2020/JASPAR2020_CORE_non-redundant_pfms_meme.txt

cd $BASEDIR/Result_NvP53/NvP53_Peaks_PseudoBulk_20220805/Explain/codetest/
tomtom -oc tomtom_t1_allr -thresh 0.1 -dist allr -no-ssc ../meme-simple.txt $MOTIFDB &>log.tomtom_conv1 &
tomtom -oc tomtom_t1_allr_ssc -thresh 0.1 -dist allr ../meme-simple.txt $MOTIFDB &>log.tomtom_conv1 &

tomtom -oc tomtom_t1_allr_trimmed -thresh 0.1 -dist allr -no-ssc ../meme-simple-trimmed.txt $MOTIFDB &>log.tomtom_conv1 &
tomtom -oc tomtom_t1_pearson_trimmed -thresh 0.1 -no-ssc ../meme-simple-trimmed.txt $MOTIFDB &>log.tomtom_conv1 &

tomtom -oc tomtom_t1_allr_trimmed_ssc -thresh 0.1 -dist allr ../meme-simple-trimmed.txt $MOTIFDB &>log.tomtom_conv1 &
tomtom -oc tomtom_t1_allr_trimmed_internal -thresh 0.1 -dist allr -no-ssc -internal ../meme-simple-trimmed.txt $MOTIFDB &>log.tomtom_conv1 &

# -------------------- annotate meme against CisTargetDB(Scenic) -------------------
MOTIFDB=$BASEDIR/Resources/scenic.meme

cd $BASEDIR/Result_NvP53/NvP53_Peaks_PseudoBulk_20220805/Explain
tomtom -oc tomtom_conv1_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme-simple.txt $MOTIFDB &>log.tomtom_conv1 &
tomtom -oc tomtom_conv1Trim_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme-simple-trimmed.txt $MOTIFDB &>log.tomtom_conv1 &
tomtom -oc tomtom_conv1Scaner_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme-scaners.txt $MOTIFDB &>log.tomtom_conv1 &
tomtom -oc tomtom_conv1ScanerTrim_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme-scaners-trimmed.txt $MOTIFDB &>log.tomtom_conv1 &

cd $BASEDIR/Result_NvP53/NvP53_Peaks_PseudoBulkT0001_20220807/Explain
tomtom -oc tomtom_conv1_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme-simple.txt $MOTIFDB &>log.tomtom_conv1 &
tomtom -oc tomtom_conv1Scaner_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme-scaners.txt $MOTIFDB &>log.tomtom_conv1 &

cd $BASEDIR/Result_NvP53/NvP53_Peaks_PseudoBulkT01_20220807/Explain
tomtom -oc tomtom_conv1_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme-simple.txt $MOTIFDB &>log.tomtom_conv1 &
tomtom -oc tomtom_conv1Scaner_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme-scaners.txt $MOTIFDB &>log.tomtom_conv1 &

# -------------------- annotate meme against CisTargetDB(Scenic) -------------------
MOTIFDB=$BASEDIR/Resources/scenic.meme

cd $BASEDIR/Result_NvP53/NvP53_AtacRna_KoWt_20220830/Explain
tomtom -oc tomtom_KoScaner_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme.KoScaner.txt $MOTIFDB &>log.meme.KoScaner.txt &
tomtom -oc tomtom_WtScaner_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme.WtScaner.txt $MOTIFDB &>log.meme.WtScaner.txt &
tomtom -oc tomtom_SharedScaner_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme.SharedScaner.txt $MOTIFDB &>log.meme.SharedScaner.txt &

cd $BASEDIR/Result_NvP53/NvP53_AtacRna_KoWt_Lung_20220830/Explain
tomtom -oc tomtom_KoScaner_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme.KoScaner.txt $MOTIFDB &>log.meme.KoScaner.txt &
tomtom -oc tomtom_WtScaner_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme.WtScaner.txt $MOTIFDB &>log.meme.WtScaner.txt &
tomtom -oc tomtom_SharedScaner_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme.SharedScaner.txt $MOTIFDB &>log.meme.SharedScaner.txt &

cd $BASEDIR/Result_NvP53/NvP53_AtacRna_KoWt_Spleen_2022902/Explain
tomtom -oc tomtom_KoScaner_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme.KoScaner.txt $MOTIFDB &>log.meme.KoScaner.txt &
tomtom -oc tomtom_WtScaner_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme.WtScaner.txt $MOTIFDB &>log.meme.WtScaner.txt &
tomtom -oc tomtom_SharedScaner_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme.SharedScaner.txt $MOTIFDB &>log.meme.SharedScaner.txt &

cd $BASEDIR/Result_NvP53/NvP53_AtacRna_KoWt_Thymus_2022902/Explain
tomtom -oc tomtom_KoScaner_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme.KoScaner.txt $MOTIFDB &>log.meme.KoScaner.txt &
tomtom -oc tomtom_WtScaner_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme.WtScaner.txt $MOTIFDB &>log.meme.WtScaner.txt &
tomtom -oc tomtom_SharedScaner_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme.SharedScaner.txt $MOTIFDB &>log.meme.SharedScaner.txt &

cd $BASEDIR/Result_NvP53/NvP53_AtacRna_KoWt_20220906/Explain
tomtom -oc tomtom_KoScaner_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme.KoScaner.axis0.txt $MOTIFDB &>log.meme.KoScaner.txt &
tomtom -oc tomtom_WtScaner_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme.WtScaner.axis0.txt $MOTIFDB &>log.meme.WtScaner.txt &
tomtom -oc tomtom_SharedScaner_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme.SharedScaner.axis0.txt $MOTIFDB &>log.meme.SharedScaner.txt &

cd $BASEDIR/Result_NvP53/NvP53_AtacRna_KoWt_FilterNum_20220907/Explain
tomtom -oc tomtom_KoScaner_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme.KoScaner.axis0.txt $MOTIFDB &>log.meme.KoScaner.txt &
tomtom -oc tomtom_WtScaner_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme.WtScaner.axis0.txt $MOTIFDB &>log.meme.WtScaner.txt &

cd $BASEDIR/Result_NvP53/NvP53_AtacRna_KoWt_P2ResNet_20220907/Explain
tomtom -oc tomtom_KoScaner_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme.KoScaner.axis0.txt $MOTIFDB &>log.meme.KoScaner.txt &
tomtom -oc tomtom_WtScaner_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme.WtScaner.axis0.txt $MOTIFDB &>log.meme.WtScaner.txt &

cd $BASEDIR/Result_NvP53/NvP53_AtacRna_KoWt_P2FilterLen_20220907/Explain
tomtom -oc tomtom_KoScaner_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme.KoScaner.axis0.txt $MOTIFDB &>log.meme.KoScaner.txt &
tomtom -oc tomtom_WtScaner_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme.WtScaner.axis0.txt $MOTIFDB &>log.meme.WtScaner.txt &
ß
cd $BASEDIR/Result_NvP53/NvP53_AtacRna_KoWt_DGE_20220908/Explain
tomtom -oc tomtom_KoScaner_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme.KoScaner.axis0.txt $MOTIFDB &>log.meme.KoScaner.txt &
tomtom -oc tomtom_WtScaner_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme.WtScaner.axis0.txt $MOTIFDB &>log.meme.WtScaner.txt &

cd $BASEDIR/Result_NvP53/NvP53_AtacRna_KoWt_P3Convs_20220908/Explain
tomtom -oc tomtom_KoScaner_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme.KoScaner.axis0.txt $MOTIFDB &>log.meme.KoScaner.txt &
tomtom -oc tomtom_WtScaner_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme.WtScaner.axis0.txt $MOTIFDB &>log.meme.WtScaner.txt &

cd $BASEDIR/Result_NvP53/NvP53_AtacRna_KoWt_P3Convs_20220908_Num16/Explain
tomtom -oc tomtom_KoScaner_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme.KoScaner.axis0.txt $MOTIFDB &>log.meme.KoScaner.txt &
tomtom -oc tomtom_WtScaner_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme.WtScaner.axis0.txt $MOTIFDB &>log.meme.WtScaner.txt &

cd $BASEDIR/Result_NvP53/NvP53_AtacRna_KoWt_P3Convs_20220914_Num4/Explain
tomtom -oc tomtom_KoScaner_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme.KoScaner.axis0.txt $MOTIFDB &>log.meme.KoScaner.txt &
tomtom -oc tomtom_WtScaner_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme.WtScaner.axis0.txt $MOTIFDB &>log.meme.WtScaner.txt &
tomtom -oc tomtom_SharedScaner_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme.SharedScaner.axis0.txt $MOTIFDB &>log.meme.SharedScaner.txt &

cd $BASEDIR/Result_NvP53/NvP53_AtacRna_KoWt_P3Convs_20220914_Num8/Explain
tomtom -oc tomtom_KoScaner_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme.KoScaner.axis0.txt $MOTIFDB &>log.meme.KoScaner.txt &
tomtom -oc tomtom_WtScaner_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme.WtScaner.axis0.txt $MOTIFDB &>log.meme.WtScaner.txt &
tomtom -oc tomtom_SharedScaner_CisTarget_t1_allr -thresh 0.1 -dist allr -no-ssc ./meme.SharedScaner.axis0.txt $MOTIFDB &>log.meme.SharedScaner.txt &

cd $BASEDIR/Result_NvP53/NvP53_AtacRna_KoWt_20220918_T0.1/Explain
tomtom -oc tomtom_KoScaner_atac -thresh 0.1 -dist allr -no-ssc meme.KoScaner.axis0.atac.txt $MOTIFDB &>log.meme.KoScaner.atac.txt &
tomtom -oc tomtom_KoScaner_rna -thresh 0.1 -dist allr -no-ssc meme.KoScaner.axis0.rna.txt $MOTIFDB &>log.meme.KoScaner.rna.txt &
tomtom -oc tomtom_WtScaner_atac -thresh 0.1 -dist allr -no-ssc meme.WtScaner.axis0.atac.txt $MOTIFDB &>log.meme.WtScaner.atac.txt &
tomtom -oc tomtom_WtScaner_rna -thresh 0.1 -dist allr -no-ssc meme.WtScaner.axis0.rna.txt $MOTIFDB &>log.meme.WtScaner.rna.txt &
tomtom -oc tomtom_SharedScaner_atac -thresh 0.1 -dist allr -no-ssc meme.SharedScaner.axis0.atac.txt $MOTIFDB &>log.meme.SharedScaner.atac.txt &
tomtom -oc tomtom_SharedScaner_rna -thresh 0.1 -dist allr -no-ssc meme.SharedScaner.axis0.rna.txt $MOTIFDB &>log.meme.SharedScaner.rna.txt &

cd $BASEDIR/Result_NvP53/NvP53_AtacRna_KoWt_P4ResNet_20220919/Explain
tomtom -oc tomtom_KoScaner_atac -thresh 0.1 -dist allr -no-ssc meme.KoScaner.axis0.atac.txt $MOTIFDB &>log.meme.KoScaner.atac.txt &
tomtom -oc tomtom_KoScaner_rna -thresh 0.1 -dist allr -no-ssc meme.KoScaner.axis0.rna.txt $MOTIFDB &>log.meme.KoScaner.rna.txt &
tomtom -oc tomtom_WtScaner_atac -thresh 0.1 -dist allr -no-ssc meme.WtScaner.axis0.atac.txt $MOTIFDB &>log.meme.WtScaner.atac.txt &
tomtom -oc tomtom_WtScaner_rna -thresh 0.1 -dist allr -no-ssc meme.WtScaner.axis0.rna.txt $MOTIFDB &>log.meme.WtScaner.rna.txt &
tomtom -oc tomtom_SharedScaner_atac -thresh 0.1 -dist allr -no-ssc meme.SharedScaner.axis0.atac.txt $MOTIFDB &>log.meme.SharedScaner.atac.txt &
tomtom -oc tomtom_SharedScaner_rna -thresh 0.1 -dist allr -no-ssc meme.SharedScaner.axis0.rna.txt $MOTIFDB &>log.meme.SharedScaner.rna.txt &

# ------------------------------------------------------------------
MOTIFDB=$BASEDIR/Resources/scenic.meme

cd $BASEDIR/Result_NvP53/NvP53_AtacRna_KoWt_Final_20221008/Explain
tomtom -oc tomtom_KoScaner_atac -thresh 0.1 -dist allr -no-ssc meme.KoScaner.axis0.atac.txt $MOTIFDB &>log.meme.KoScaner.atac.txt &
tomtom -oc tomtom_KoScaner_rna -thresh 0.1 -dist allr -no-ssc meme.KoScaner.axis0.rna.txt $MOTIFDB &>log.meme.KoScaner.rna.txt &
tomtom -oc tomtom_WtScaner_atac -thresh 0.1 -dist allr -no-ssc meme.WtScaner.axis0.atac.txt $MOTIFDB &>log.meme.WtScaner.atac.txt &
tomtom -oc tomtom_WtScaner_rna -thresh 0.1 -dist allr -no-ssc meme.WtScaner.axis0.rna.txt $MOTIFDB &>log.meme.WtScaner.rna.txt &
tomtom -oc tomtom_SharedScaner_atac -thresh 0.1 -dist allr -no-ssc meme.SharedScaner.axis0.atac.txt $MOTIFDB &>log.meme.SharedScaner.atac.txt &
tomtom -oc tomtom_SharedScaner_rna -thresh 0.1 -dist allr -no-ssc meme.SharedScaner.axis0.rna.txt $MOTIFDB &>log.meme.SharedScaner.rna.txt &
