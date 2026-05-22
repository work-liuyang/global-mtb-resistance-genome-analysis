export PATH=/mnt/project/software/mamba_env_liuyang/envs/circos/bin:$PATH
#export PATH=/work/share/acgzcbkvwv/Software/conda_envs/svgtools/bin:$PATH

for i in Australia China India Moldova USA
do
  [ ! -d "$i/prepare_data" ] && mkdir -p "$i/prepare_data"
  cd $i && python3 ../make_karyotype.py ../${i}_decade_snp.xls && \
  circos -conf ../general.conf -outputdir ./ -outputfile ${i}_circos && \
  /mnt/project/software/mamba_env_liuyang/bin/rsvg-convert -f pdf -o ${i}_circos.pdf ${i}_circos.svg 
  cd ../
done
