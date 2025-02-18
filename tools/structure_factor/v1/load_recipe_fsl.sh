#!/bin/bash

USER=$(echo $HOME | cut -d'/' -f3);

#added by Anaconda3 5.3.0 installer
# >>> conda init >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$(CONDA_REPORT_ERRORS=false
'/fslhome/${USER}/fsl_groups/fslg_softmatter/bin/anaconda3/bin/conda'
shell.bash hook 2> /dev/null)"
if [ $? -eq 0 ]; then
    \eval "$__conda_setup"
else
    if [ -f
"/fslhome/${USER}/fsl_groups/fslg_softmatter/bin/anaconda3/bin/etc/profile.d/conda.sh"
]; then
        .
"/fslhome/${USER}/fsl_groups/fslg_softmatter/bin/anaconda3/bin/etc/profile.d/conda.sh"
        CONDA_CHANGEPS1=false conda activate base
    else
        \export
PATH="/fslhome/${USER}/fsl_groups/fslg_softmatter/bin/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda init <<<
