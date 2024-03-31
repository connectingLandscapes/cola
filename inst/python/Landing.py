import streamlit as st
import base64

# Set page config
st.set_page_config(
    page_title="COLA",
    page_icon="🌍",
    layout="centered",
    initial_sidebar_state="collapsed"
)

with open("./style.css") as f:
    st.markdown(f'<style>{f.read()}</style>', unsafe_allow_html=True)
    
# Page title
st.title("COLA")

# Quote
st.markdown(
    """
    > "There is no Planet B 🌍."
    """,
)
st.sidebar.title("About")
st.sidebar.info(
    """
    - Web App URL: <https://rbachati-cdpop-webapp-app-ugs378.streamlit.app/>
    - GitHub repository: <https://github.com/rbachati/cdpop_webapp>
    """
)

st.sidebar.title("Contact")
st.sidebar.info(
    """
    Patrict Jantz
    [GitHub](https://github.com/forest-rev) |[LinkedIn]()
    """
)


# Additional Text Content
st.markdown(
    """
    ## What is COLA?
    COLA is a web application dedicated to promoting awareness about environmental sustainability.
    We provide data-driven insights and resources to help individuals and organizations make more eco-friendly decisions.
    """
)

def add_bg_from_local(image_file):
    with open(image_file, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read())
    st.markdown(
    f"""
    <style>
    .stApp {{
        background-image: url(data:image/{"jpg"};base64,{encoded_string.decode()});
        background-size: cover
    }}
    </style>
    """,
    unsafe_allow_html=True
    )
add_bg_from_local('./images/background.jpg')

st.image("./images/background.png", use_column_width=True)

row1_col1, row1_col2 = st.columns(2)
with row1_col1:
    st.image("./images/kavango.gif")
    st.image("./images/200w.webp", use_column_width=True)

with row1_col2:
    st.image("https://github.com/giswqs/data/raw/main/timelapse/goes.gif")
    st.image("https://github.com/giswqs/data/raw/main/timelapse/fire.gif")
    

# More links or navigation
st.markdown(
    """
    ---
    ##### Explore More
    - [AUthor/Team Section](#)
    - [Our Mission](#)
    - [Contact Us](#)
    - [Source Repository](https://github.com/rbachati/cdpop_webapp)
    """,
)