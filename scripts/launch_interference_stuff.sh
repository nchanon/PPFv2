#./bin/histograms_creator m_dilep 2017
echo "  1) rootfile creation for combine"
python ./bin/combine_unrolled.py $1 $2
echo "  2) datacard creation for combine"
./bin/card_creator $1 $2 Interference
echo "  3) export to lyoserv"
python ./scripts/export_combine.py interference $2