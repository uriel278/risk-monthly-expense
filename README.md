# Will I be able to cover my monthly expenses?

The project in this repository explores models to quantify the risk of a person not being able to cover their monthly expenses. The data source is the 2021 National Survey of Financial Inclusion, carried out by INEGI (Instituto Nacional de Estadistica y Geografia). 

The main objectives of this work are (i) the identification of covariates that can explain the risk of interest, and (ii) the development
of a model that can be used to predict the risk based on new observations of the set of covariates. The covariates included in the model are of different nature such as binary, integer, categorical, ordinal and continous. These covariates aim to capture different categories such as financial awareness and socieconomic proxies.

The modeling strategy consists of a regularized logistic regression model. The regularization penalty corresponds to the elastic net (L1 + L2). With this approach, the model discards the non relevant covariates by setting them to zero. Three different link functions for the logistic regression models were compared.

Predictive performance was assessed by AUC plots for each of the three link functions. The results indicate similar performance and the logit link, which provides with easier interpretability, was chosen for the final model.
