copy 180314 into ~/privatemodules/GITM/
module load GITM/180314
cd /short/n23/<path to your run directory>
makeRunDir.sh 
cp $GITM_HOME/rungitm.pbs ./
qsub rungitm.pbs
