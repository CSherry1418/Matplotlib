#!/usr/bin/env python
# coding: utf-8

# ## Observations and Insights 

# 

# In[23]:


# Dependencies and Setup
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st
import numpy as np

# Study data files
mouse_metadata_path = "mouse_metadata.csv"
study_results_path = "mouse_study_results.csv"

# Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path)

# Combine the data into a single dataset
combined_df = pd.merge(mouse_metadata, study_results)

# Display the data table for preview
combined_df = pd.merge(mouse_metadata, study_results, on="Mouse ID")
combined_df.head()


# In[24]:


# Checking the number of mice.
mouse_metadata["Mouse ID"].value_counts()
study_results["Mouse ID"].value_counts()


# In[25]:


# Getting the duplicate mice by ID number that shows up for Mouse ID and Timepoint. 
mouse_duplicates = combined_df.loc[combined_df.duplicated(subset=['Mouse ID', 'Timepoint',]),'Mouse ID'].unique()
mouse_duplicates


# In[26]:


# Create a clean DataFrame by dropping the duplicate mouse by its ID.
combined_df_rem = pd.DataFrame.drop_duplicates(combined_df)
combined_df_rem.head()


# In[27]:


# Checking the number of mice in the clean DataFrame.
combined_df_rem["Mouse ID"].value_counts()


# ## Summary Statistics

# In[28]:


# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen
# Use groupby and summary statistical methods to calculate the following properties of each drug regimen: 
# mean, median, variance, standard deviation, and SEM of the tumor volume. 
# Assemble the resulting series into a single summary dataframe.

mean = combined_df_rem.groupby('Drug Regimen')['Tumor Volume (mm3)'].mean()
median = combined_df_rem.groupby('Drug Regimen')['Tumor Volume (mm3)'].median()
variance = combined_df_rem.groupby('Drug Regimen')['Tumor Volume (mm3)'].var()
standard_deviation = combined_df_rem.groupby('Drug Regimen')['Tumor Volume (mm3)'].std()
sem = combined_df_rem.groupby('Drug Regimen')['Tumor Volume (mm3)'].sem()
summary_df = pd.DataFrame({"Mean": mean, "Median": median, "Variance": variance, "Standard Deviation": standard_deviation, "SEM": sem})
summary_df
 


# In[29]:


# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen
# Using the aggregation method, produce the same summary statistics in a single line


# ## Bar and Pie Charts

# In[30]:


# Generate a bar plot showing the total number of measurements taken on each drug regimen using pandas.
grouped_df = pd.DataFrame(combined_df.groupby(["Drug Regimen"]).count()).reset_index()
regimen_measurements_pan = grouped_df[["Drug Regimen", "Mouse ID"]]
regimen_measurements_pan = regimen_measurements_pan.rename(columns={"Mouse ID": "Count"})
regimen_measurements_pan = regimen_measurements_pan.set_index("Drug Regimen")
regimen_measurements_pan.plot(kind="bar", figsize=(10,3))
plt.title("Total Number per Drug Regimen")

plt.show()


# In[66]:


# Generate a bar plot showing the total number of measurements taken on each drug regimen using pyplot.

regimen_measurements_pyp = summary_df.index.tolist()
regimen_count = (combined_df.groupby(["Drug Regimen"])["Age_months"].count()).tolist()
x_axis = regimen_measurements_pyp
x_axis = np.arange(len(regimen_count))
plt.figure(figsize=(10,3))
plt.bar(x_axis, regimen_count)
plt.title("Total Number per Drug Regimen")
plt.xlabel("Drug Regimen")
tick_locations = [value for value in x_axis]
plt.xticks(tick_locations, regimen_measurements_pyp, rotation=90)


# In[74]:


# Generate a pie plot showing the distribution of female versus male mice using pandas
gender_df = pd.DataFrame(combined_df.groupby(["Sex"]).count()).reset_index()
gender_df = gender_df[["Sex", "Mouse ID"]]
gender_df = gender_df.rename(columns={"Mouse ID": "Count"})
plt.figure(figsize=(10,6))
ax1 = plt.subplot(121, aspect='equal')
gender_df.plot(kind='pie', y = "Count", ax=ax1, autopct='%1.1f%%', 
 startangle=90, shadow=False, labels=gender_df['Sex'], legend = False)


# In[73]:


# Generate a pie plot showing the distribution of female versus male mice using pyplot
gender_count = (combined_df.groupby(["Sex"])["Age_months"].count()).tolist()
labels = ["Male", "Female"]
plt.pie(gender_count, labels=labels, autopct="%1.1f%%", shadow=False, startangle=90)


# ## Quartiles, Outliers and Boxplots

# In[77]:


# Calculate the final tumor volume of each mouse across four of the treatment regimens:  
# Capomulin, Ramicane, Infubinol, and Ceftamin

# Start by getting the last (greatest) timepoint for each mouse
final_tumor_vol_df = combined_df.groupby(['Drug Regimen', 'Mouse ID']).last()['Timepoint']

# Merge this group df with the original dataframe to get the tumor volume at the last timepoint
drug_reg_merge_df = pd.merge(final_tumor_vol_df, combined_df, on=("Mouse ID", "Timepoint"),how="left")
drug_reg_merge_df


# In[14]:


# Put treatments into a list for for loop (and later for plot labels)


# Create empty list to fill with tumor vol data (for plotting)


# Calculate the IQR and quantitatively determine if there are any potential outliers. 

    
    # Locate the rows which contain mice on each drug and get the tumor volumes
    
    
    # add subset 
    
    
    # Determine outliers using upper and lower bounds
    


# In[15]:


# Generate a box plot of the final tumor volume of each mouse across four regimens of interest


# ## Line and Scatter Plots

# In[93]:


# Generate a line plot of tumor volume vs. time point for a mouse treated with Capomulin
capomulin_df = combined_df.loc[combined_df["Drug Regimen"] == "Capomulin"]
capomulin_df = capomulin_df.reset_index()
capomulin_mouse = capomulin_df.loc[capomulin_df["Mouse ID"] == "s185"]
capomulin_mouse = capomulin_mouse.loc[:, ["Timepoint", "Tumor Volume (mm3)"]]
capomulin_mouse = capomulin_mouse.reset_index(drop=True)
capomulin_mouse.set_index('Timepoint').plot(figsize=(10,10))


# In[95]:


# Generate a scatter plot of average tumor volume vs. mouse weight for the Capomulin regimen
average_cap = pd.DataFrame(capomulin_mouse.groupby(["Mouse ID", "Weight (g)"])["Tumor Volume (mm3)"].mean()).reset_index()
average_cap = average_cap.rename(columns={"Tumor Volume (mm3)": "Average Volume"})
average_cap = average_cap.set_index('Mouse ID')
average_cap.plot(kind="scatter", x="Weight (g)", y="Average Volume", grid=True, figsize=(4,4),
              title="Weight Vs. Average Tumor Volume")
plt.show()


# ## Correlation and Regression

# In[18]:


# Calculate the correlation coefficient and linear regression model 
# for mouse weight and average tumor volume for the Capomulin regimen


# In[ ]:




