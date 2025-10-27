#!/bin/bash
# Christian Panse <cp@fgcz.ethz.ch>
#
#SBATCH --mem-per-cpu=4096M
#SBATCH --cpus-per-task=4
#SBATCH -N 1 # count of nodes                   
#SBATCH -n 1 # cpu count / number of tasks
#SBATCH --partition=fcc
#SBATCH --nodelist=fgcz-c-072
#SBATCH --job-name=rawrr
#SBATCH -o /tmp/slurm.%x.%j.out
#SBATCH -e /tmp/slurm.%x.%j.err

set -e
set -x

input=$1
output=$2


now=$(date "+%s") 
inputdate=$(date -r $input "+%s")

#[[ $(($now - $inputdate)) -lt 36000 ]]
    

[[ -f $output ]] && { echo "Warning: $output exists"; exit 0; }

mkdir -p $(dirname $output)
echo $input $output

R -q --no-save <<EOF
.autoQC <- function(f){
  iRTmz <- c(487.2571, 547.2984, 622.8539, 636.8695, 644.8230, 669.8384, 683.8282,
    683.8541, 699.3388, 726.8361, 776.9301)

  names(iRTmz) <- c('LGGNEQVTR', 'YILAGVENSK', 'GTFIIDPGGVIR', 'GTFIIDPAAVIR',
    'GAGSSEPVTGLDAK', 'TPVISGGPYEYR', 'VEATFGVDESNAK',
    'TPVITGAPYEYR', 'DGLDAASYYAPVR', 'ADVTPADFSEWSK',
    'LFLQFGAQGSPFLK')

    rawrr::readChromatogram(f, mass = iRTmz, tol = 10, type = 'xic', filter = 'ms') -> C
    C |> plot(diagnostic = TRUE)
}
fi <- '$input'
fo <- '/tmp/slurm.rawrr.$SLURM_JOB_ID.pdf'
paste0('rawrr version', packageVersion('rawrr')) -> st
pdf(fo, 19,12)
rawrr::readChromatogram(fi, type = 'bpc') |> plot()
rawrr::readChromatogram(fi, type = 'tic') |> plot()
if (grepl('autoQC', fi)){
  .autoQC(fi)
}

fi |> rawDiag::readRaw() -> RAW
RAW |> rawDiag::plotChargeState()
RAW |> rawDiag::plotCycleTime()
RAW |> rawDiag::plotLockMassCorrection()
RAW |>rawDiag::plotMzDistribution()
RAW |>rawDiag::plotScanTime()
RAW |>rawDiag::plotCycleLoad()
RAW |>rawDiag::plotInjectionTime()
RAW |>rawDiag::plotMassDistribution()
RAW |>rawDiag::plotPrecursorHeatmap()
RAW |>rawDiag::plotTicBasepeak()

plot(0, 0, type='n', xlab='', ylab=''); text(0,0,st, cex=3)
dev.off()
EOF

## scp $SSHLOGIN:/tmp/orawrr.pdf r35:$output
test -s /tmp/slurm.rawrr.$SLURM_JOB_ID.pdf && cp -v /tmp/slurm.rawrr.$SLURM_JOB_ID.pdf $output

