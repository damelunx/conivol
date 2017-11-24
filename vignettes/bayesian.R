## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "conic-intrinsic-volumes_figures/"
)

## ----load-pkgs, include=FALSE--------------------------------------------
library(conivol)
library(tidyverse)
library(knitr)
library(png)

## ----load-imgs, include=FALSE--------------------------------------------
img_paths <- list( dl="bayes_diagrams/bayes_direct_logconc.png",
                   dnl="bayes_diagrams/bayes_direct_nologconc.png",
                   idle="bayes_diagrams/bayes_indirect_logconc_even.png",
                   idlo="bayes_diagrams/bayes_indirect_logconc_odd.png",
                   idnle="bayes_diagrams/bayes_indirect_nologconc_even.png",
                   idnlo="bayes_diagrams/bayes_indirect_nologconc_odd.png" )
img_dpi <- 450

## ----disp-dnl, echo=FALSE------------------------------------------------
include_graphics(img_paths$dnl, dpi=img_dpi)

## ----disp-dl, echo=FALSE-------------------------------------------------
include_graphics(img_paths$dl, dpi=img_dpi)

## ----disp-idnle, echo=FALSE----------------------------------------------
include_graphics(img_paths$idnle, dpi=img_dpi)

## ----disp-idnlo, echo=FALSE----------------------------------------------
include_graphics(img_paths$idnlo, dpi=img_dpi)

## ----disp-idle, echo=FALSE-----------------------------------------------
include_graphics(img_paths$idle, dpi=img_dpi)

## ----disp-idlo, echo=FALSE-----------------------------------------------
include_graphics(img_paths$idlo, dpi=img_dpi)

