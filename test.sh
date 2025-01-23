#before running the test pipeline, please decompress all the .gz files from "./test_data/"
perl ./IFRAMLOG_v0.1.pl  --vcf ./test_data/Query_test.vcf --ref ./test_data/Query_test_ref.fa --out_dir ./out/ --out Test.flk150 --OG_ref ./test_data/Outgroup_test_ref.fa
