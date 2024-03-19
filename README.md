# AMPActiPred: a three-stage computational framework for predicting antibacterial peptides and their activity levels

[**AMPActiPred**](https://awi.cuhk.edu.cn/~AMPActiPred) is a three-stage prediction framework based on a deep forest architecture designed to  **identify ABP and characterize their targets across bacterial strains, culminating in the prediction of ABP activity levels** . Given the potential functional activity of ABPs against one or more different bacteria strains, we establish a one-against-all multitask classification strategy to predict peptides with multiple targets. We employed the deep forest algorithm to construct the model, leveraging a hierarchical architecture similar to deep neural networks, capable of effectively processing and mining original features for enhanced predictive performance. Additionally, AMPActiPred incorporates several efficient peptide descriptors that effectively represent the composition and physicochemical properties of peptides. Furthermore, we conducted a feature interpretability analysis to assist researchers in identifying crucial features influencing antibacterial activity.  **AMPActiPred is the first computational framework capable of predicting targets and activity levels of ABPs** , which demonstrates convincing performance in characterizing ABPs and their targets while accurately predicting ABP activity levels, offering a promising avenue for expediting the development of novel ABPs.

---

## Requirements

AMPActiPred has been successfully tested in a `Python 3.6` environment and relies on the following packages. Please ensure that the following packages are installed.

```
biopython==1.76
modlamp==4.3.0
pandas==1.1.5
numpy==1.18.5
deep-forest==0.1.5
joblib==1.2.0
```

## Usage Guide

First, download the models for three stages from the following links respectively:

* First stage: [https://awi.cuhk.edu.cn/~AMPActiPred/download/first_stage_model_final.tar.gz](https://awi.cuhk.edu.cn/~AMPActiPred/download/first_stage_model_final.tar.gz)
* Second stage: [https://awi.cuhk.edu.cn/~AMPActiPred/download/second_stage.tar.gz](https://awi.cuhk.edu.cn/~AMPActiPred/download/second_stage.tar.gz)
* Third stage: [https://awi.cuhk.edu.cn/~AMPActiPred/download/third_stage.tar.gz](https://awi.cuhk.edu.cn/~AMPActiPred/download/third_stage.tar.gz)

Then execute the following command to extract:

```js
tar zxvf first_stage_model_final.tar.gz
tar zxvf second_stage.tar.gz
tar zxvf third_stage.tar.gz
```

Then execute predict.py to run the prediction program. For example:

```
python predict.py -i sequenceText.fasta -o record.json --Ecoli True --Saureus True --Ecoli True -predSpecie True --predActivity True
```

The input of the program is a Fasta file containing stored protein sequences. To facilitate further processing by the computer, the output of the model is a file in JSON format. For more information about parser arguments, please refer to the `predict.py`. Below are explanations for some key parser arguments.

```js
-i The path of input file
-o The path of outputfile
-predSpecie Select whether to predict functional activity against specific bacteria
--predActivity Select whether to predict the activity level of ABP
```

## Other auxiliary files and scripts

* The dataset used for the model can be found here: [https://github.com/lantianyao/AMPActiPred/tree/main/data](https://github.com/lantianyao/AMPActiPred/tree/main/data)
* [extract_feature.ipynb](https://github.com/lantianyao/AMPActiPred/blob/main/extract_feature.ipynb) provides code to **extract peptide descriptors** used in this study.
* [use the trained model.ipynb](https://github.com/lantianyao/AMPActiPred/blob/main/use%20the%20trained%20model.ipynb) is a quick tutorial telling you how to  **use a trained deep forest model** .

For details of deep forest python package, please refer to the API [here](https://deep-forest.readthedocs.io/en/latest/ "here").

If you have any questions, please feel free to contact me.

**Name: Lantian Yao Email: [lantianyao@link.cuhk.edu.cn](mailto:lantianyao@link.cuhk.edu.cn)**

