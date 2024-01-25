#! /bin/bash

# Make movie from initial frames of simulation file.
# Requires FFmpeg.

# DEFAULT PARAMETERS

__PYTHON_DEFAULT=$(which python)    # default python executable
__FFMPEG_DEFAULT=$(which ffmpeg)    # default FFmpeg executable

# HELP

usage() {
cat <<< "
Make movie from initial frames of simulation file using module \`cells.read\`.

Requires FFmpeg. (see https://ffmpeg.org/download.html)

SYNOPSIS

    [bash] $(basename "${BASH_SOURCE[0]}") [OPTIONS] SIMULATION_FILE_NAME

OPTIONS

    -h  Display this help.

    -r  Make rainbow plot.

    -p  Python executable.
        DEFAULT: $__PYTHON_DEFAULT
    -f  FFmpeg executable.
        DEFAULT: $__FFMPEG_DEFAULT

    -y  \"Yes\" to all FFmpeg prompts.
"
}

# OPTIONS

while getopts "hrp:f:y" OPTION; do
    case $OPTION in
        r)  # rainbow plot
            __RAINBOW=true;;
        p)  # python executable
            __PYTHON=$OPTARG;;
        f)  # FFmpeg executable
            __FFMPEG=$OPTARG;;
        y)  # yes to FFmpeg
            __YES=true;;
    esac
done
shift $(expr $OPTIND - 1)

if [[ -z "$@" ]]; then
    usage;
    exit 1;
fi

__PYTHON=${__PYTHON-$__PYTHON_DEFAULT}
__FFMPEG=${__FFMPEG-$__FFMPEG_DEFAULT}

# SCRIPT

tmp_dir=$(mktemp -d)

cat >> ${tmp_dir}/movie.py <<EOF
"""
Save frames from simulation file.
"""

from cells.read import Read, _progressbar

import os

# simulation file
r = Read("$(realpath $1)")
# plotted frames
frames = list(r.t0)

_progressbar(0)
for frame in frames:

    # plot
    r.plot(frame ${__RAINBOW:+, rainbow=r.t0[0]})
    # save
    r.fig.savefig(os.path.join("$tmp_dir", "%05d.png" % frames.index(frame)))

    _progressbar((frames.index(frame) + 1)/len(frames))

EOF

# save frames
$__PYTHON ${tmp_dir}/movie.py
# make movie
$__FFMPEG ${__YES:+-y} -r 5 -f image2 -s 1280x960 \
    -i ${tmp_dir}/%5d.png -pix_fmt yuv420p ${1%.*}.mp4
# make movie with H.265 compression (not compatible with all players)
# requires FFmpeg configured with --enable-gpl --enable-libx265
# $__FFMPEG ${__YES:+-y} -r 5 -f image2 -s 1280x960 -pix_fmt yuv420p \
#     -i ${tmp_dir}/%5d.png -vcodec libx265 -crf 28 ${1%.*}.mp4

