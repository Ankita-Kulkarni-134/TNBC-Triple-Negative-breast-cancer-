import os
import requests
from playwright.sync_api import sync_playwright
import pandas as pd
import time
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("admetlab_processing.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


def process_file(file_path, page, output_dir, max_retries=3):
    retries = 0
    while retries < max_retries:
        try:
            logger.info(f"Processing file: {file_path} (Attempt {retries + 1}/{max_retries})")

            # Get the base filename without extension and add .csv
            base_filename = os.path.basename(file_path)
            base_name_without_ext = os.path.splitext(base_filename)[0]
            output_filename = f"{base_name_without_ext}.csv"
            download_path = os.path.join(output_dir, output_filename)

            # Navigate to the page with increased timeout
            logger.info("Navigating to ADMETlab screening page...")
            page.goto('https://admetlab3.scbdd.com/server/screening', timeout=300000)
            page.wait_for_load_state('networkidle', timeout=300000)

            # Locate the file input element
            file_input = page.wait_for_selector('input[type="file"]', timeout=60000)
            logger.info("File input found.")
            file_input.set_input_files(file_path)
            logger.info(f"File '{file_path}' set in the input element.")

            # Locate and click the submit button
            submit_button = page.wait_for_selector('button.btn.btn-success[type="button"]', timeout=60000)
            logger.info("Submit button found.")
            submit_button.scroll_into_view_if_needed()

            # Use JavaScript click for more reliability
            page.evaluate('(button) => button.click()', submit_button)

            # Wait for the results to be processed
            logger.info("Waiting for results to be processed...")
            page.wait_for_timeout(90000)  # Increased to 90 seconds to be safer

            # Get the current page URL after processing
            result_url = page.url
            logger.info(f"Result URL: {result_url}")

            # Look for the "Download as CSV" button in the top-right corner
            logger.info("Looking for 'Download as CSV' button...")
            # Try multiple selector strategies
            download_button = None
            selectors = [
                'button:has-text("Download as CSV")',
                'a:has-text("Download as CSV")',
                '.btn-success:has-text("CSV")',
                'button.btn-download',
                '.dropdown-menu a:has-text("CSV")'
            ]

            for selector in selectors:
                try:
                    download_button = page.wait_for_selector(selector, timeout=5000)
                    if download_button:
                        logger.info(f"Found download button with selector: {selector}")
                        break
                except:
                    continue

            if not download_button:
                logger.error("Could not find the download button with standard selectors")
                logger.info("Taking screenshot to debug...")
                page.screenshot(path=f"debug_screenshot_{int(time.time())}.png")

                # Last resort: try to find all buttons on the page
                buttons = page.query_selector_all('button, a.btn')
                logger.info(f"Found {len(buttons)} buttons/links on the page")
                for i, btn in enumerate(buttons):
                    text = btn.inner_text()
                    logger.info(f"Button {i}: {text}")
                    if "csv" in text.lower() or "download" in text.lower():
                        download_button = btn
                        logger.info(f"Found likely download button with text: {text}")
                        break

            if not download_button:
                raise Exception("Could not find download button")

            # Start waiting for download
            with page.expect_download(timeout=60000) as download_info:
                # Click the download button
                download_button.click()
                logger.info("Clicked download button")

            # Wait for download to complete and save the file
            download = download_info.value
            logger.info(f"Download started: {download.suggested_filename}")

            # Save to the specified path with the desired filename
            download.save_as(download_path)
            logger.info(f"File saved to {download_path}")

            # Record the successful download
            success_log = os.path.join(output_dir, "successful_downloads.txt")
            with open(success_log, 'a') as f:
                f.write(f"{file_path}\t{download_path}\t{result_url}\n")

            # Successfully processed, break out of retry loop
            return download_path

        except Exception as e:
            logger.error(f"Error: {str(e)}")
            retries += 1
            if retries < max_retries:
                logger.info(f"Retrying in 20 seconds... ({retries}/{max_retries})")
                time.sleep(20)  # Longer wait between retries
            else:
                logger.error(f"Maximum retries reached for {file_path}. Moving to next file.")
                # Record the failure
                failed_log = os.path.join(output_dir, "failed_downloads.txt")
                with open(failed_log, 'a') as f:
                    f.write(f"{file_path}\n")
                return None


def run(playwright):
    # Disable SSL warnings for requests
    try:
        import urllib3
        urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
    except:
        pass

    # Launch browser with SSL errors ignored
    browser = playwright.chromium.launch(
        headless=False,
        args=['--ignore-certificate-errors']
    )

    # Create context with SSL errors ignored
    context = browser.new_context(
        ignore_https_errors=True,
        accept_downloads=True  # Enable downloads
    )
    page = context.new_page()

    # Output directory
    output_dir = r"C:/Users/adina/OneDrive/Desktop/admet/Rushi/sample_admet/sample_admet_out"
    os.makedirs(output_dir, exist_ok=True)

    # Files to process
    file_paths = [


r"C:/Users/adina/OneDrive/Desktop/admet/Rushi/split_outs/Combined_unique_1500.txt",
r"C:/Users/adina/OneDrive/Desktop/admet/Rushi/split_outs/Combined_unique_1501.txt",






    ]

    # Process each file
    results = []
    for file_path in file_paths:
        download_path = process_file(file_path, page, output_dir)
        if download_path:
            results.append((file_path, download_path))
            logger.info(f"Successfully processed {file_path}")
        else:
            logger.warning(f"Failed to process {file_path}")

    # Summary
    logger.info(f"Processing complete. Successful: {len(results)}, Failed: {len(file_paths) - len(results)}")

    browser.close()


if __name__ == "__main__":
    with sync_playwright() as playwright:
        run(playwright)