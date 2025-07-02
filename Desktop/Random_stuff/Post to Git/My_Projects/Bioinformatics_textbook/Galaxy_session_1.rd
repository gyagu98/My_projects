# -- Getting Started with Galaxy 

# Genrietta Yagudayeva

# 2025-03-21

## Session 33: Getting Started with Galaxy

## 33.1 Introduction to Galaxy

[**Galaxy**](https://usegalaxy.org/){target="_blank"} is a *free* and open-source web platform that enables users to analyze biological data. [**Galaxy**](https://usegalaxy.org/){target="_blank"} is designed to be a graphical user interface (GUI)-friendly, it is a web-based platform made for bioinformatics analysis. It is most known for enabling users to perform computational biology tasks without requiring programming skills.

Here are some key features of [**Galaxy**](https://usegalaxy.org/){target="_blank"}:
1) Allows programming style analysis without command-line tools
GUI for ease of use
2) Supports various bioinformatics applications, including biostatistics, NGS data processing, RNA-Seq analysis, and variant calling
3) Facilitates data analysis, workflow authoring, tool integration, training and education, and infrastructure management
4) Enables users to create and share automated workflows to generate easier scientific reproducibility

Let's explore the Galaxy Architecture!
```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/1_galaxy_architecture.png")
```

There are three panels to take note of: The *left* panel shows a menu where you can upload data, available tools, pre-made workflows, histories, and pages. The *middle* panel is the primary space where users will execute their analysis and work. The *right* panel is the history panel, which will track and organize your data analysis and allow users to revisit and modify workflows.

In summary, Galaxy is an excellent option for better management of in silico analysis, especially for novices in the computational sciences.

## 33.2 Hands-on Activity: Let’s navigate Galaxy!

Follow this Hands-on Activity to navigate in Galaxy:
1) Step 1: Make a Galaxy account. To use Galaxy, you must make an account. In the upper right corner of the homepage, select Login or Register > Register here. Choose an email you have access to
After registering, you will be prompted to verify your account via email. 

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/2_galaxy_login.png")
```
After routine verification, you should have free access to Galaxy with 250 GB of free storage (you can always upgrade for more)

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/3_galaxy_login.png")
```

2) Step 2: Explore all the menu options

*Upload*

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/4_galaxy_upload.png")
```

This will be where you upload the raw data that you can either choose from your device, from remote options (bioinformatics databases), or by physically pasting data. It will take some time to upload your data, depending on the size of the file.

*Tools*

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/5_galaxy_tools.png")
```

Tools are software packages, often command-line applications, that are wrapped to allow users to interact with them through a user-friendly, form-driven interface, enabling reproducible analysis. As of today, March 10, 2025, the Galaxy ToolShed hosts over 9500 distinct and modular software packages (tools)
You can group and filter them by function, tasks, scientific discipline, workflows, etc. There are general text tools, genomic file manipulation tools, standard genomics tools, genomics analysis tools, statistics and visualization tools, genomic toolkits, and domain tools. 

*Workflows*

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/6_workflows.png")
```

Workflows are sequences of tools and dataset actions that run in a defined order, enabling complex analyses to be performed automatically and reproducibly. Galaxy Workflows facilitate rigorous, reproducible analysis pipelines that can be shared with the community. You can either make a workflow, share them, or use publicly available ones.

*Workflow Invocations*

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/7_workflow_invocations.png")
```

Workflow Invocations are a single run or execution of a defined workflow, which is a series of tools and dataset actions executed in sequence.

*Visualizations*

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/8_visualizations.png")
```

Visualizations are tools and methods that allow users to visually represent and analyze data, particularly genomic data, using various chart types, genome browsers, and other specialized views.


*Histories*

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/9_histories.png")
```

Histories are chronological records of your analyses, containing all datasets uploaded or generated, along with the tools and parameters used, enabling reproducibility and tracing of data provenance. Galaxy uses histories to organize your datasets and analysis steps, acting as a workspace where you can upload, run tools, and store results. These are vital as you work through projects; whether in collaboration or juggling multiple projects, it is essential that before running any analysis, you first name a history.

*Histories Multiview*

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/10_histories_multiview.png")
```

Histories Multiview serves as an organization tool to view and manipulate multiple histories simultaneously, allowing for easier comparison and organization of analyses by enabling the dragging of datasets between them.

*Datasets*

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/11_datasets.png")
```

Datasets are the inputs and outputs of each step in an analysis project in Galaxy. They serve as structured collections of data, typically representing input, intermediate, or output files used or produced during an analysis, which can be uploaded, manipulated, and analyzed using Galaxy's tools.

*Pages*

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/12_pages.png")
```

Pages are the primary way to communicate your Galaxy analyses so that other people can easily view, reproduce, or extend your analyses.

*Libraries*

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/13_libraries.png")
```

Libraries are a way to store and share datasets within a group of Galaxy users or with everyone who has access to a specific Galaxy instance.

*Notifications*

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/14_notifications.png")
```

Right now, you won’t have any notifications, but notifications are handy to notify you when certain file uploading is done, data analysis is completed, or if you work with other collaborators, you may benefit from UseGalaxy communications as you exchange information/data/workflows and more.

So please, go ahead, make an account and explore on your own! Using the many tools and menu options of Galaxy, you can craft projects on our own and allow for a more friendly approach to computational analyses

## Session 33.3 What is Galaxy Training, and why should we all do it?

Let’s go ahead and explore [**Galaxy Training**](https://training.galaxyproject.org/){target="_blank"}, which is a collection of tutorials developed and maintained by the worldwide Galaxy community. 

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/15_training.png")
```

Select Learning Pathways from the top toolbar, or visit this [**link**](https://training.galaxyproject.org/training-material/learning-pathways/){target="_blank"}, you will see sets of tutorials curated for you by community experts to form a coherent set of lessons around a topic, building up knowledge as you go.

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/16_learning_pathways.png")
```

Galaxy training designers recommend to follow the tutorials in the order they are listed in the pathway. For our bioinformatics focus, let’s focus specifically on the training designed for scientists. An example learning pathway that is suggested for all users to start with is *Introduction to Galaxy and Sequence Analysis*, can be found on this [**link**](https://training.galaxyproject.org/training-material/learning-pathways/intro-to-galaxy-and-genomics.html){target="_blank"}.

Each learning pathway will explain to you why you would need to learn it, which modules to learn, how long each module would take, and the difficulty level. So, in the case of the  *Introduction to Galaxy and Sequence Analysis*, you would learn:
The basics of Galaxy and how to use Galaxy for analysis (Module 1 of approximately 1 hour and 40 minutes). The most common first steps of any genome analysis are quality control and mapping or assembly of your genomic sequences (Module 2 of approximately 5 hours).

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/17_learning_pathways_2.png")
```

You can choose from many learning pathways specifically for your project by doing the following: in the upper right corner, click “Search Tutorials” specific to your interest. For example, if I am studying poison dart frogs and want to determine how their genes are differentially expressed and how they differ between species and/or organs, I would type “transcriptomics” into the search bar.

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/18_learning_pathways_3.png")
```

Now, I have several resources to read through and find which one is most applicable to me. You specifically want to find the ones that say “Hands-on” (to learn how to process your data from scratch and possibly design your workflow), “Workflow” (to process your data through an already published workflow), or “Slides” (to learn the theoretical background about your experiment). 

Galaxy training is vital for all people using bioinformatics for  a few different reasons:
1) You will learn bioinformatics more quickly as it provides a platform to learn and practice various bioinformatics concepts and tools, with a wide variety of tutorials and resources available
2) You will use skills specific to your project by following specific learning pathways, modules, and tutorials
3) You will collect a good grasp on the ability of Galaxy and its tools to unlock its potential fully

You should review Galaxy training in your own time to determine what pathway(s) you want to take for your project.

## 33.4 Hands-on Activity: Uploading files, understanding metadata, and running a simple workflow

In this Hands-on Activity, you will learn how to create/use histories on Galaxy, upload data, and run a simple workflow, so please follow these steps:

*Step 1 Creating a New History and Uploading Data in Galaxy*:
a) Create a new history by navigating to the History Panel on the right side of the Galaxy interface
b) Click on the plus sign icon and select Create New History

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/19_hands-on_1a.png")
```

c) Rename the history to “my_first_galaxy_run” by clicking on the pencil, entering the new name, and hitting Save

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/20_hands-on_1b.png")
```

d) Download this sample FASTQ file on this link [**here**](https://zenodo.org/record/3977236/files/female_oral2.fastq-4143.gz){target="_blank"} titled “female_oral2.fastq-4143.gz” (This is a microbiome sample from a snake Jacques et al. 2021 taken from the QC Galaxy Tutorial written by Cameron Hyde et al. 2016; see tutorial here on this link [**here**](https://training.galaxyproject.org/training-material/topics/sequence-analysis/tutorials/quality-control/tutorial.html){target="_blank"})

A FASTQ file is a text-based format used in bioinformatics to store both a biological sequence (like DNA or RNA) and its corresponding quality scores, essentially providing information about the confidence level of each base call made during sequencing, making it the standard format for raw data from high-throughput sequencing machines like Illumina platforms, the file you downloaded is a tiny file size sample of what a FASTQ file can look like. Make sure to unzip the file and open the FASTQ file via a text editor software, like Sublime or BBedit. 

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/21_hands-on_1c.png")
```

The basic anatomy of a FASTQ file is the following: 
Each read consists of four lines:
Line 1: Sequence Identifier, which starts with @
Line 2: Raw DNA Sequence, which are the nucleotide bases (A, T, C, G)
Line 3: Separator Line, which is always just a + symbol
Line 4: Quality Scores, which are ASCII-encoded quality values
Higher ASCII values means Higher sequencing quality (less sequencing error probability)

Go back to Galaxy, click on the Upload Data button (top left panel), and ensure you are in your proper history “my_first_galaxy_run”.

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/22_hands-on-1d.png")
```

Typically, you could choose one of the following methods to upload files (as mentioned previously), but in this case, we will select the local file you just downloaded, “female_oral2.fastq-4143.gz”. Set the file format to FASTQ (if not auto-detected). 

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/23_hands-on-1e.png")
```

Click Start to begin the upload. Wait for the file to appear in the history panel (right side of the screen) and ensure it turns green (indicating successful upload). 

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/24_hands-on-1f.png")
```

You should now see the file in your History on the right panel

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/25_hands-on-1g.png")
```

*Step 2 Understanding Metadata in Galaxy to familiarize yourself with dataset metadata and learn how to modify it for clarity and proper workflow execution*: 

a) Press the uploaded file in the history to see the immediate metadata

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/26_hands-on-2a.png")
```

The immediate metadata tells us:
Dataset Number and Name:
7 indicates that this is the 7th dataset in my history
"female_oral2.fastq-4143" is the name of your uploaded file

Dataset Status (Green Background): 
The green color means the dataset has been successfully uploaded and processed without errors. If the dataset were still processing, it would appear in yellow. If an error had occurred, it would turn red.

Add Tags:
Clicking this allows you to assign tags to your dataset for easier identification in workflows. This is useful when juggling analyses and/or projects.

Dataset Information:
812 sequences indicates that this FASTQ file contains 812 sequencing reads, format fastqsanger confirms the detected file format is FASTQ. database with the question mark means the dataset is not yet associated with a reference genome or database

Dataset Description Field:
The text field says uploaded fastq file. This is a user-defined description. You can edit this by clicking inside the box and adding notes.

Action Icons Below the Description Box performs a specific function:
- Save Icon means to *Download* the dataset to your local computer.

- Link Icon means *Share or Link* Information and generate a shareable link or view dataset details.

- i Icon means *Information* and gives a total metadata background on the file and provides a command line view of the code written (for learning programming directly).

- Bar graph Icon means *View Dataset Statistics* to check summary statistics, file size, and sequence details.

- Folder Icon means *Workflow Connections* to view dataset relationships within workflows and analyses.

Now, on the top right of the file name, the action icons says:

- Eye Icon to *To view the dataset* and clicking this opens a detailed view of the dataset in Galaxy (it should be the same as what you saw when you opened it on your local computer)

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/27_hands-on-2b.png")
```

- Pencil Icon means to *Edit attributes* edit dataset attributes like name and info, add an annotation, or add a database/build

- Trash Icon means to *Delete Dataset* and deletes the dataset from your history. It moves the dataset to the deleted items section (it can still be recovered if needed)

*Step 3 Running a Simple Workflow in Galaxy by executing a basic bioinformatics workflow:* 

Part 3A: Quality Control tools

1) Running FastQC (Quality Control Check).
2) Locate the FastQC Tool.
3) In the left panel, search for FastQC under the tools section.
4) Click on FastQC: Read Quality Control.

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/28_hands-on3a.png")
```

5) Select “Raw read data from your current history” and choose the “sample_reads.fastq” file. Note here: Galaxy will tell you which file formats are accepted for the specified tool. 

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/29_hands-on3b.png")
```

6) On the selected tool, you will be able to see (with examples from the FASTQC tool):

List of various parameters you can select/deselect depending on your goal (for instance, contaminant list, adapter list, etc.) – for usual practice, keep the settings at default

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/30_hands-on3c.png")
```

Additional options (for example, email notification)

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/31_hands-on3d.png")
```

Help section to explain why the tool is used, about the tool, inputs/outputs, tutorials, and a help forum).

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/32_hands-on3e.png")
```

7) Now, Click “Run Tool” to run the analysis.

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/33_hands-on3f.png")
```

Note that your History will add the job(s).

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/34_hands-on3g.png")
```


Wait for the job to complete (the dataset is orange during the run, green when complete, and red if there is an error).

8) Analyze both FastQC Outputs

*FastQC on data Webpage* Click the eye icon to visualize your webpage data.

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/35_hands-on.png")
```

*FastQC on data RawData* Click the eye icon to visualize your raw data.

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/36_hands-on.png")
```

9) Review key metrics in the Webpage file: 

*Per-base sequence quality*
Shows the quality scores of nucleotides at each position in the reads.High-quality reads have scores above 20-30 (Phred scale)
A downward trend at the end of the reads suggests degradation and may require trimming.

*GC content*
Represents the proportion of Guanine (G) and Cytosine (C) bases.
It should be uniform for a given species or sequencing run.
GC bias (unexpected peaks or drops) can indicate contamination or sequencing errors.

*Adapter contamination*
Checks for adapter sequences left in the reads. Adapters need to be removed using Trimmomatic or Cutadapt before further analysis
High levels of adapter contamination may indicate short reads or over-sequencing

*Identify any potential quality issues*
Low per-base quality means to consider trimming poor-quality bases, GC-content anomalies means possible contamination or sequencing bias, high adapter content means to need trimming before alignment or downstream processing.

Part 3B: Trimming tools
A good next step after QC is trimming your data to clean it up before further analysis.

So let’s use Cutadapt to clean up our data and see if our summary metrics improve.

Steps:
1) Locate the Cutadapt Tool
2) Search for Cutadapt in the tools panel
3) Click on Cutadapt: Remove adapter sequences

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/37_hands-on.png")
```

4) In Cutadapt, select: Single-ended reads, input your original file, in Read 1 Adapter, in 3’ End Adapters, select custom sequence, and in the following box, type: CTGTCTC

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/38_handson.png")
```

In Other Read Trimming Options select Quality cutoff(s) (R1): 22, and in Read Filtering Option,s select Minimum length (R1): 15

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/39_handson.png")
```

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/40_handson.png")
```

5) Hit the Run Tool and let the file run, and wait for the file to turn green

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/41_handson.png")
```

6) Now, use this output file to see if the FASTQC report changed by plugging this file into the FASTQC tool and hitting Run Tool

7) Let’s hit the eye icon on our data file of the FASTQC web page and observe

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/42_handson.png")
```

Based on the previous original data, our summary for FASTQC is better as we now do not have the adapter, and the base quality is better and improved.


Step 3C: Putting it all together and running a Simple Workflow in Galaxy by executing a basic bioinformatics workflow


As you can imagine, bioinformaticians deal with a lot of data simultaneously. Wouldn’t it be nice to automate our pipelines and have data processed in bulk? That’s where workflows come in, and we can utilize our tools from Galaxy in an automated fashion! Quality control and trimming sequences are also typical tasks a bioinformatician would take in cleaning their data before further analysis. So, let’s design a simple FASTQC + Cutadapt pipeline, by following these steps:

1) Navigate to Workflows and click the Create button

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/43_workflow_handson.png")
```

2) Label it “my_first_galaxy_workflow

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/44_workflow.png")
```


3) Select the tools you used from the previous analyses (FASTQC and Cutadapt) with all the same parameters and input them into your new Galaxy workflow

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/45_workflow.png")
```


4) Click the save icon at the top next to your title, and now you have your first workflow!

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/46_workflow.png")
```


5) Click “Run,” and now you can run your very own workflow without the need to search for any tools; you can also modify it to run the steps, too

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/47_workflow.png")
```

6) Click “Run Workflow” to see your jobs getting done simultaneously and result in the same as your previous steps

```{r echo = FALSE, out.width = "100%"}
knitr::include_graphics("/Users/genriettayagudayeva/Desktop/Galaxy_session_1/Session_1/Images/48_workflow.png")
```

You can repeat your workflow as often as needed to perfect your analyses

Congratulations, you finished the first session of Galaxy Learning! In the next session, we will explore building a specific project.
