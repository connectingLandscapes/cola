# streamlit_app.py
import streamlit_authenticator
import streamlit as st
from ipysheet import from_dataframe, to_dataframe
import ipysheet
import streamlit as st
import pandas as pd
import os
import glob

# Set the title of the app
st.set_page_config(page_title="CDPOP Simulation Model", page_icon=":bar_chart:", layout="wide")
st.title("LandScape Genetics (CDPOP)")
    
if streamlit_authenticator.check_password():
    # Define function to run CDPOP script
    def run_cdpop(input_file_path, output_folder):
        # Change directory to where the CDPOP script is located
        os.chdir(os.path.join(os.path.dirname(__file__), "CDPOP", "src"))
        
        # Only keep the file name (strip away the directory)
        input_file_name = os.path.basename(input_file_path)

        # Run the CDPOP script with input and output paths
        os.system(f"python CDPOP.py ../data {input_file_name} ../data/{output_folder}")

    def main():       
        # Define input and output file paths
        input_file = st.sidebar.file_uploader("Upload input CSV file", type=["csv"])
        if input_file is not None:
        # Save the uploaded file to disk
            input_path = os.path.join(os.path.dirname(__file__), "CDPOP", "data", input_file.name)
            st.session_state.input_df = pd.read_csv(input_file)
            
        output_folder = st.sidebar.text_input("Enter output folder name. Ex: CDPOP_Scenario1")

        # Define the editable columns
        editable_cols = [
            "xyfilename",
            "output_years",
            "output_unicor",
            "dispcdmat",
            "matemovethresh",
            "Fdispmovethresh",
            "Mdispmovethresh",
        ]
        
        # Define new column names
        new_col_names = {
            "xyfilename": "XY Locations file",
            "output_years": "generations/years of saved genotypes",
            "output_unicor": "xy locations for output years",
            "dispcdmat": "dispersal cost distance matrix",
            "matemovethresh": "mating distance threshold",
            "Fdispmovethresh": "Female Dispersion Threshold",
            "Mdispmovethresh": "Male Dispersion Threshold",
        }

        # If input file is specified, load and edit the input CSV file
        if 'input_df' in st.session_state:

            # Add a checkbox to control the visibility of the table
            show_df = st.checkbox("select the checkbox for the parameters file")

            # If the checkbox is checked, display the table
            if show_df:
                st.write(st.session_state.input_df)

            st.header("Edit input parameter file")
            # Make a copy of your dataframe with only the columns you want to be editable
            editable_df = st.session_state.input_df[editable_cols].copy()
            
            #Rename the columns
            editable_df.rename(columns=new_col_names, inplace=True)

            # Use the experimental data editor
            edited_df = st.data_editor(editable_df, num_rows="dynamic", key="data_editor")

            # Handle any changes made in the data editor
            if st.session_state.data_editor:
                # Convert the new column names back to the original names
                edited_df.rename(columns={v: k for k, v in new_col_names.items()}, inplace=True)
                # Update the original dataframe with the changes made in the data editor
                st.session_state.input_df.update(edited_df)

            if st.button("Save changes"):
                st.session_state.input_df.to_csv(input_path, index=False)
                st.success("Changes saved successfully!")

            # Add a Run button and a progress bar
            if st.button("Run"):
                with st.spinner('Running CDPOP...'):
                    run_cdpop(input_path, output_folder)
                st.success(f"Done running CDPOP. The output was saved in the folder {output_folder} in the data directory.")
                
            if st.sidebar.button("List Output Files"):
                # Parent directory
                data_dir = os.path.join(os.path.dirname(__file__), "CDPOP", "data")

                # Find all matching directories
                matching_dirs = glob.glob(os.path.join(data_dir, output_folder + "*"))

                if matching_dirs:  # if list is not empty
                    # Collect all file names
                    if 'all_files' not in st.session_state:
                        st.session_state.all_files = []
                        for dir in matching_dirs:
                            files = os.listdir(dir)
                            for file in files:
                                st.session_state.all_files.append(f"{dir}: {file}")

                                # Display the files in a dropdown in the sidebar
                                selected_file = st.sidebar.selectbox("Select an Output File:", st.session_state.all_files)
                                st.sidebar.write(f"You selected: {selected_file}")
                            else:
                                st.sidebar.write(f"No output directories starting with {output_folder} found.")


    if __name__ == "__main__":
        main()