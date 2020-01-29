“ODD for the IBM `Behaviour-Community.R`”
================
“Team KeithChaets”
29 January, 2020

Here we outline the ODD protocol for the IBM `Behaviour-Community.R`,
based on Grimm et al 2006 & 2010. The model described …

-----

## 1\. Overview

-----

  - **Purpose**
      - The model simulates aggressive interactions between reef fishes,
        and how the nature of these interactions influence community
        structure \[29/01/19\]
      - The model tracks individual fish (agents) from multiple species.
        State variables that describe each agent are:
          - species \_ baseline probability of aggression, consumption,
            movement, reproduction and mortality
          - age \_ accrue through model time steps
          - sex \_ Assigned through Bernoulli process
          - energetic scope \_ Informed by bioenergetics models, and
            emerges as a result of behavioural processes at each time
            step
  - Process overview and scheduling

-----

## 2\. Design concepts

-----

  - design concepts

-----

## 3\. Details

-----

  - Initialization
  - Input
  - Submodels
