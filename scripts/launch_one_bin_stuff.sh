#./bin/histograms_creator m_dilep 2017
echo "  1) rootfile creation for combine"
python ./bin/combine_one_bin.py $1 $2
echo "  2) datacard creation for combine"
./bin/card_creator $1 $2 OneBin
echo "  3) export to lyoserv"
python ./scripts/export_combine.py one_bin $2