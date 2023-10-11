# MICB475_Group13
### authors: Juliet Malkowski, Lina Anwari, Abigail Cho, Màiri MacAulay, Betty Hong

# Project Proposal
### Proposed title: 
-often reflects the overarching goal of your project

### Introduction and Background: 
-provides the premise for why/how your dataset was collected in addition to an overview of the studies that have already been conducted on your dataset or relevant to your dataset (consider what you presented for your P2 review oral presentation)

### Research Objectives: 
-explains the areas of interest that you wish to pursue including specific research questions and corresponding hypotheses

### Experimental Aims and Rationale: 
-list of aims that your team hopes to achieve in order to address your research objectives and briefly describes how each aim will help achieve the research objectives

### Proposed Approach: 
-tabular summary of the purpose and proposed approach for each experimental aim

| AIM 1A: Does sheet washing frequency impact skin microbial composition? | AIM 1B: Does shower recency impact skin microbial composition? | AIM 2: Do gender differences in hygiene practices impact skin microbial composition? | AIM 3: Does time spent outside impact skin microbial composition in conjunction with hygiene practices?|
|---|---|---|---|
|1A-1: Edit metadata document to divide sheet washing frequency into high/low|1B-1: Edit metadata document to divide last shower day into recent/not recent|2-1: Import and demultiplex raw sequence reads|3-1: Edit metadata document to divide hours outside into high/low 
|1A-2: Import and demultiplex raw sequence reads (demux.qza)|1B-2: Import and demultiplex raw sequence reads|2-2: Perform sequence quality control and generate a feature table (table.qza)|3-2: Import and demultiplex raw sequence reads|
|1A-3: Perform sequence quality control and generate a feature table (table.qza)|1B-3: Perform sequence quality control and generate a feature table (table.qza)|2-3: Filter metadata by gender|3-3: Perform sequence quality control and generate a feature table (table.qza)|
|1A-4: Filter metadata for only samples that contain a value for sheet washing frequency|1B-4: Filter metadata for only samples that contain a value for shower recency|2-4: Filter metadata for samples that contain a value for either sheet washing frequency or shower recency|3-4: Filter metadata for samples that have a value for time outside|
|1A-5: Generate a filtered table|1B-5: Generate a filtered table|2-5: Generate a filtered table|3-5: Filter metadata for samples that contain a value for either sheet washing frequency or shower recency|
|1A-6: Run alpha and beta diversity metrics and analyze for interesting results|1B-6: Run alpha and beta diversity metrics and analyze for interesting results|2-6: Run alpha and beta diversity metrics and analyze for interesting results|3-6: Generate a filtered table|1A-7: Run taxonomic analysis to determine strain differences between groups (generate taxonomic bar blot)|1B-7: Run taxonomic analysis to determine strain differences between groups (generate taxonomic bar blot)|2-7: Run taxonomic analysis to determine strain differences between groups (generate taxonomic bar blot)|3-7: Run alpha and beta diversity metrics and analyze for interesting results|
|1A-7: Run taxonomic analysis to determine strain differences between groups (generate taxonomic bar blot)|1B-7: Run taxonomic analysis to determine strain differences between groups (generate taxonomic bar blot)|2-7: Run taxonomic analysis to determine strain differences between groups (generate taxonomic bar blot)|3-7: Run alpha and beta diversity metrics and analyze for interesting results|
|1A-8: Filter to remove unwanted ASVs|1B-8: Filter to remove unwanted ASVs|2-8: Filter to remove unwanted ASVs|3-8: Run taxonomic analysis to determine strain differences between groups (generate taxonomic bar blot)|
|1A-9: Export data and perform statistical analysis on alpha and/or beta diversity metrics and/or differential abundance testing for taxonomic metrics using R|1B-9: Export data and perform statistical analysis on alpha and/or beta diversity metrics and/or differential abundance testing for taxonomic metrics using R|2-9: Export data and perform statistical analysis on alpha and/or beta diversity metrics and/or differential abundance testing for taxonomic metrics using R|3-9: Filter to remove unwanted ASVs|
||||3-10: Export data and perform statistical analysis on alpha and/or beta diversity metrics and/or differential abundance testing for taxonomic metrics using R|

- I'll put in the names of the actual columns/groups we make
- and I can make what alpha/beta diversity metrics we choose more specific

### Overview Flowchart: 
-visual representation of research objectives or questions, corresponding experimental aims, and corresponding analysis/approach.

### Weekly Timeframe:
```mermaid
gantt
dateFormat  YYYY-MM-DD
title Group 13 Project 2 Timeline

section Week 6[Oct.9-Oct.13]
Update Project Proposal on Github                               :done,          first_1,    2023-10-02, 2023-10-10
Write title,introduction,research objectives,experimental aims  :active,        first_2,    2023-10-08, 2023-10-10
Weekly Question Agenda on Proposed Approach/Overview            :active,        first_3,    2023-10-08, 2023-10-10
Meeting with Evelyn/Chris                                       :milestone,     m1,         2023-10-11,1d

section Week 7[Oct.16-Oct.20]
Do processing steps in QIIME2 and flowchart overview            :active,        first_4,    2023-10-11, 2023-10-17
Weekly Question Agenda on QIIME2 processing                     :active,        first_5,    2023-10-16, 2023-10-18
Meeting with Evelyn/Chris                                       :milestone,     initial,    2023-10-18,1d
Project Proposal due                                            :milestone,     initial,    2023-10-22,1d

section Week 8[Oct.23-Oct.27]
Weekly Question Agenda on                                       :active,        first_6,    2023-10-19, 2023-10-25
Meeting with Evelyn/Chris                                       :milestone,     initial,    2023-10-25,1d

section Week 9[Oct.30-Nov.3]
Weekly Question Agenda on                                       :active,        first_7,    2023-10-26, 2023-11-01
Meeting with Evelyn/Chris                                       :milestone,     initial,    2023-11-01,1d


section Week 10[Nov.6-Nov.10]
Weekly Question Agenda on                                       :active,        first_8,    2023-11-02, 2023-11-08
Meeting with Evelyn/Chris                                       :milestone,     initial,    2023-11-08,1d


section Week 11[Nov.13-Nov.17]
Midterm Break                                                   :done,          break,    2023-11-13, 2023-11-17

section Week 12[Nov.20-Nov.24]
Weekly Question Agenda on                                       :active,        first_9,    2023-11-14, 2023-11-22
Meeting with Evelyn/Chris                                       :milestone,     initial,    2023-11-22,1d


section Week 13[Nov.27-Dec.1]
Weekly Question Agenda on                                       :active,        first_10,    2023-11-23, 2023-11-29
Meeting with Evelyn/Chris                                       :milestone,     initial,    2023-11-29,1d


section Week 14[Dec.4-Dec.8]
Oral Presentations                                              :milestone,    initial,    2023-12-04,1d
Draft Manuscript due                                            :milestone,    initial,    2023-12-10,1d


section Week 15[Dec.11-Dec.15]
Final Manuscript & Lab Notebook due                             :milestone,    initial,    2023-12-17,1d

```

### Dataset Overview: 
-In order to complete this section, you will need to complete the processing steps in QIIME2 (up until the rarefaction curve) to extract the following information and describe your dataset. Please use this checklist Download checklistto ensure you describe these elements within your proposal and include the summary table and 2 figures listed in the checklist. 

### Participation Report:
-A breakdown of each team member’s contribution to preparing the proposal
Juliet: coded weekly timeframe

### References:
-follow ASM referencing guidelinesLinks to an external site.

1- mermaid r diagram to use for weekly timeline: https://mermaid.js.org/ 
