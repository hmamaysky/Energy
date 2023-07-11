#!/bin/bash

cd /shared/share_mamaysky-glasserman/energy_drivers/2023
rm -rf DataProcessing_acc
rm -rf DataProcessing_mod
mkdir -p DataProcessing_acc/article_measure
mkdir -p DataProcessing_mod/article_measure

cd DataProcessing_acc/article_measure
ln -s ../../DataProcessing/article_measure/3gram 3gram
ln -s ../../DataProcessing/article_measure/4gram 4gram
ln -s ../../DataProcessing/article_measure/entropy entropy
ln -s ../../DataProcessing/article_measure/sentiment sentiment
ln -s ../../DataProcessing/article_measure/total total
mkdir topic_allocation

cd ../../DataProcessing_mod/article_measure
ln -s ../../DataProcessing/article_measure/3gram 3gram
ln -s ../../DataProcessing/article_measure/4gram 4gram
ln -s ../../DataProcessing/article_measure/entropy entropy
ln -s ../../DataProcessing/article_measure/sentiment sentiment
ln -s ../../DataProcessing/article_measure/total total
mkdir topic_allocation