echo "##########################"
echo "## PyPlotFramework v2.0 ##"
echo "##########################"

sleep 1
echo " > Results directory creation"
mkdir -p results/{2016,2017,2018}/{comparaison,flattree,stats,other,systematics}
mkdir -p results/pheno/base

sleep 1
echo " > Inputs directory creation"
mkdir -p inputs/{2016,2017,2018}/{DATA,MC,JEC,ALT}
mkdir -p inputs/{pheno,timed,other}

sleep 1
echo " > Combine directory creation"
mkdir -p combine/{2016,2017,2018}/{sme,interference,one_bin,unrolled}/{inputs,results}

sleep 1
echo "Installtion finished ! :)"