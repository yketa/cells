#! /bin/bash

# Make movie from initial frames of simulation file.
# Requires FFmpeg.

# DEFAULT PARAMETERS

__PYTHON_DEFAULT=$(which python)    # default python executable
__FFMPEG_DEFAULT=$(which ffmpeg)    # default FFmpeg executable

# HELP

usage() {
cat <<< "
Make movie from initial frames of simulation file using modules \`cells.read\`
and \`cells.plot\`.

Requires FFmpeg. (see https://ffmpeg.org/download.html)

SYNOPSIS

    [bash] $(basename "${BASH_SOURCE[0]}") [OPTIONS] SIMULATION_FILE_NAME
    [bash] $(basename "${BASH_SOURCE[0]}") [OPTIONS] -d FRAMES_DIRECTORY_NAME

OPTIONS

    -h  Display this help.

    -r  Make rainbow plot.
        NOTE: Does nothing if -d option is used.
    -c  Clear the plot of all cell colouring.
        NOTE: Does nothing if -d option is used.
              Options -r and -c are exclusionary.
    -f  Display forces on vertices.

    -d  Make movie from frames in directory.
        DEFAULT: (not specified)
        NOTE: Frames must be *.png and ordered.

    -p  Python executable.
        DEFAULT: $__PYTHON_DEFAULT
    -F  FFmpeg executable.
        DEFAULT: $__FFMPEG_DEFAULT

    -y  \"Yes\" to all FFmpeg prompts.
"
}

# OPTIONS

while getopts "hrcfd:p:F:y" OPTION; do
    case $OPTION in
        h)  # help
            usage; exit 1;;
        r)  # rainbow plot
            __RAINBOW=true;;
        c)  # clear plot
            __CLEAR=true;;
        f)  # force plot
            __FORCES=true;;
        d)  # frames directory
            __DIR="$OPTARG";;
        p)  # python executable
            if [ "$OPTARG" != "" ]; then __PYTHON=$OPTARG; fi;;
        F)  # FFmpeg executable
            if [ "$OPTARG" != "" ]; then __FFMPEG=$OPTARG; fi;;
        y)  # yes to FFmpeg
            __YES=true;;
    esac
done
shift $(expr $OPTIND - 1)

if [[ -z "$1" && -z "$__DIR" ]]; then
    usage;
    exit 1;
fi

__PYTHON=${__PYTHON-$__PYTHON_DEFAULT}
__FFMPEG=${__FFMPEG-$__FFMPEG_DEFAULT}

# SCRIPT

# make frames
if [[ -z "$__DIR" ]]; then

    __DIR=$(mktemp -d);

    cat >> ${__DIR}/movie.py <<EOF
"""
Save frames from simulation file.
"""

from cells.read import Read, _progressbar
${__FORCES:+from cells.plot import plot_forces}

import os

# simulation file
r = Read("$(realpath $1)")
# plotted frames
frames = list(r.t0)

_progressbar(0)
for frame in frames:

    # plot
    r.plot(frame
        ${__RAINBOW:+, rainbow=r.t0[0]} ${__FORCES:+, override=plot_forces})
    # save
    r.fig.savefig(os.path.join("$__DIR", "%05d.png" % frames.index(frame)))

    _progressbar((frames.index(frame) + 1)/len(frames))

EOF

    # save frames
    $__PYTHON ${__DIR}/movie.py;

fi

# make movie
__MOVIE=${1:-movie}
__MOVIE=${__MOVIE%.*}.mp4
$__FFMPEG ${__YES:+-y} -r 5 -f image2 -s 1280x960 \
    -pattern_type glob -i "${__DIR}/*.png" -pix_fmt yuv420p $__MOVIE
# make movie with H.265 compression (not compatible with all players)
# requires FFmpeg configured with --enable-gpl --enable-libx265
# $__FFMPEG ${__YES:+-y} -r 5 -f image2 -s 1280x960 -pix_fmt yuv420p \
#     -pattern_type glob -i "${__DIR}/*.png" -vcodec libx265 -crf 28 $__MOVIE

