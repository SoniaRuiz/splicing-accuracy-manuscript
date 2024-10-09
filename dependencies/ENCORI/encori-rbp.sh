#!/bin/bash
# Loop each RBP and download corresponding iCLIP data from the ENCORI platform (https://rnasysu.com/encori/)

names='A1CF AARS AATF ABCF1 ABT1 ACO1 ADAR ADD1 AGGF1 AGO1 AIMP1 AKAP1 AKAP8 AKAP8L ALYREF APOBEC3C AQR ASCC1 ATP5C1 ATXN1 ATXN2 AUH BCCIP BCLAF1 BOLL BOP1 BUD13 CALR CASC3 CCAR1 CCAR2 CCDC124 CCDC86 CCNL1 CDC40 CEBPZ CELF1 CIRBP CKAP4 CMTR1 CNOT4 CNOT7 CORO1A CPEB1 CPEB4 CPSF6 CPSF7 CSTF1 CSTF2 CSTF2T DAZ3 DAZAP1 DDX1 DDX19B DDX20 DDX21 DDX24 DDX27 DDX28 DDX3X DDX3Y DDX42 DDX47 DDX5 DDX51 DDX52 DDX55 DDX59 DDX6 DEAF1 DGCR8 DHX29 DHX30 DHX33 DKC1 DNAJC17 DNAJC2 DNAJC21 DROSHA EEF2 EFTUD2 EIF2S1 EIF2S2 EIF3A EIF3D EIF3G EIF3H EIF4A2 EIF4A3 EIF4B EIF4G1 EIF4G2 EIF4H ELAVL4 ESF1 ESRP1 ETF1 EWSR1 EXOSC5 EXOSC9 FAM120A FASTKD1 FASTKD2 FIP1L1 FKBP4 FMR1 FTO FUBP1 FUBP3 FUS FXR1 FXR2 G3BP1 G3BP2 GEMIN5 GLRX3 GNL3 GOLGB1 GPKOW GRSF1 GRWD1 GTF2F1 HDGF HLTF HMGB1 HNRNPA0 HNRNPA1 HNRNPA2B1 HNRNPAB HNRNPC HNRNPCL1 HNRNPD HNRNPDL HNRNPF HNRNPH1 HNRNPH2 HNRNPK HNRNPL HNRNPLL HNRNPM HNRNPU HNRNPUL1 HSPD1 IGF2BP1 IGF2BP2 IGF2BP3 ILF2 ILF3 KHDRBS1 KHDRBS2 KHDRBS3 KHSRP KIF1C KRR1 LARP1 LARP4 LARP7 LIN28B LSM11 LSM4 MAGOH MAK16 MARK2 MATR3 MBNL1 METAP2 MSI1 MSI2 MTPAP NAA15 NCBP2 NELFE NFX1 NIP7 NKRF NOL12 NOLC1 NONO NOVA1 NPM1 NSUN2 NUFIP2 NUP35 NUPL2 NUSAP1 NXF1 PA2G4 PABPC1 PABPC3 PABPC4 PABPN1 PABPN1L PARN PCBP1 PCBP2 PCBP3 PCBP4 PES1 PHF6 PKM PNN PNPT1 POLK POLR2G PPIG PPIL4 PPP1R8 PRPF4 PRPF6 PRPF8 PRR3 PRRC2C PSIP1 PSPC1 PTBP1 PTBP3 PUF60 PUM1 PUM2 PUS1 QKI RACK1 RALY RAVER1 RBFOX2 RBFOX3 RBM14 RBM15 RBM15B RBM17 RBM22 RBM23 RBM24 RBM25 RBM27 RBM34 RBM39 RBM4 RBM41 RBM45 RBM47 RBM4B RBM5 RBM6 RBMS2 RBMS3 RBMX2 RC3H1 RCC2 RECQL RECQL5 RPL23A RPLP0 RPS10 RPS11 RPS19 RPS2 RPS24 RPS3 RPS3A RPS5 RRP9 RTF1 SAFB2 SART3 SBDS SCAF4 SERBP1 SF1 SF3A3 SF3B1 SF3B4 SFPQ SLBP SLC4A1AP SLTM SMN1 SMNDC1 SND1 SNRNP200 SNRNP70 SNRPA SRFBP1 SRP68 SRPK2 SRSF1 SRSF10 SRSF11 SRSF2 SRSF3 SRSF4 SRSF5 SRSF7 SRSF8 SRSF9 SSB SSRP1 STAU1 STIP1 SUB1 SUCLG1 SUGP2 SUPT6H SUPV3L1 SYNCRIP TAF15 TARDBP TBRG4 TFIP11 TIA1 TIAL1 TNRC6A TOE1 TRA2A TRIM56 TRIP6 TRNAU1AP TROVE2 TUFM U2AF1 U2AF2 UBAP2L UBE2L3 UCHL5 UNK UPF1 UPF2 UTP18 UTP3 WDR3 WDR43 WRN XPO1 XPO5 XRCC5 XRCC6 XRN1 XRN2 YBX3 YTHDC2 YWHAG ZC3H11A ZC3H8 ZCRB1 ZFC3H1 ZFP36 ZNF106 ZNF326 ZNF622 ZRANB2 SAFB AGO2 NIPBL ZNF800 SDAD1 STAU2'
celltypes='HepG2 K562'

for name in $names
do

  for celltype in $celltypes
  do
  
    echo $celltype
    curl "https://rna.sysu.edu.cn/encori/api/RBPTarget/?assembly=hg38&geneType=mRNA&RBP=$name&clipExpNum=1&pancancerNum=0&target=all&cellType=$celltype" > "ENCORI_hg38_RBPTarget_"$name"_"$celltype".txt"
   
  done

done