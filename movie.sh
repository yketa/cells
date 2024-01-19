#! /bin/bash

# Make movie from initial frames of simulation file.
# Requires ffmpeg.

while getopts "rp:f:y" OPTION; do
    case $OPTION in
        r)  # rainbow plot
            __RAINBOW=true;;
        p)  # python executable
            __PYTHON=$OPTARG;;
        f)  # ffmpeg executable
            __FFMPEG=$OPTARG;;
        y)  # yes to ffmpeg
            __YES=true;;
    esac
done
shift $(expr $OPTIND - 1)

__PYTHON=${__PYTHON-python} # python executable
__FFMPEG=${__FFMPEG-ffmpeg} # ffmpeg executable

tmp_dir=$(mktemp -d)

cat >> ${tmp_dir}/movie.py <<EOF
"""
Save frames from simulation file.
"""

from cells.vm import Read, _progressbar

import os

r = Read("$(realpath $1)")                                                      # simulation file
frames = list(r.t0)                                                             # plotted frames
for frame in frames:
    r.plot(frame ${__RAINBOW:+, rainbow=r.t0[0]})                               # plot
    r.fig.savefig(os.path.join("$tmp_dir", "%05d.png" % frames.index(frame)))   # save
    _progressbar((frames.index(frame) + 1)/len(frames))

EOF

$__PYTHON ${tmp_dir}/movie.py
$__FFMPEG ${__YES:+-y} -r 5 -f image2 -s 1280x960 -pix_fmt yuv420p \
    -i ${tmp_dir}/%5d.png -vcodec libx265 -crf 28 ${1%.*}.mp4

