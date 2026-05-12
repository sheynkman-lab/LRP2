/*
 * FRAGPIPE_AUTHENTICATE: Interactive license acceptance and tool download
 *
 * This module handles the one-time setup for FragPipe with THREE modes:
 *
 * MODE 1 - INTERACTIVE (Default, like FragNFlow):
 *   - Checks if tools already exist
 *   - If not: displays license, collects user info, sends registration, prompts for token
 *   - If yes: skips download, uses cached tools
 *
 * MODE 2 - TOKEN-ONLY (Simplified for HPC):
 *   - User provides only --fragpipe_token parameter
 *   - Skips registration (user already ran curl command to get token)
 *   - Goes straight to downloading tools with provided token
 *
 * MODE 3 - FULL NON-INTERACTIVE (Legacy HPC mode):
 *   - User provides all parameters via command line
 *   - No prompts, fully automated
 *
 * Based on FragNFlow implementation: https://github.com/ronalabrcns/FragNFlow
 */

//***************
//****COLORS*****
//***************
RED = "\u001B[31m"
GREEN = "\u001B[32m"
YELLOW = "\u001B[33m"
CYAN = "\u001B[36m"
RESET = "\u001B[0m"
//***************

process FRAGPIPE_AUTHENTICATE {
    tag "FragPipe setup"
    label 'process_single'
    storeDir "${params.fragpipe_tools_dir}"

    input:
    val first_name      // Can be null for interactive mode
    val last_name       // Can be null for interactive mode
    val email           // Can be null for interactive mode
    val institution     // Can be null for interactive mode
    val token           // Can be null for interactive mode
    val license_accept  // Can be false for interactive mode

    output:
    path "msfragger/MSFragger*.jar", emit: msfragger_jar
    path "ionquant/IonQuant*.jar", emit: ionquant_jar
    path "diatracer/diaTracer*.jar", emit: diatracer_jar
    val true, emit: ready
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Detect mode:
    // - MODE 1 (Interactive): All params empty/false
    // - MODE 2 (Token-only): Only token provided, rest empty/false
    // - MODE 3 (Full non-interactive): All params provided
    def has_token = token != ''
    def has_user_info = first_name != '' && last_name != '' && email != '' && institution != ''
    def has_license = license_accept == true

    def interactive_mode = !has_token && !has_user_info && !has_license
    def token_only_mode = has_token && !has_user_info && !has_license
    def full_mode = has_token && has_user_info && has_license
    """
    #!/bin/bash
    set -euo pipefail

    # Terminal colors
    RED="\\033[31m"
    GREEN="\\033[32m"
    YELLOW="\\033[33m"
    CYAN="\\033[36m"
    RESET="\\033[0m"

    echo ""
    echo "\${CYAN}=========================================================================="
    echo "                    FRAGPIPE AUTHENTICATION & SETUP"
    echo "==========================================================================\${RESET}"
    echo ""

    mkdir -p msfragger ionquant diatracer

    #
    # STEP 1: Check if tools already exist
    #
    TOOLS_EXIST=false
    if ls msfragger/MSFragger*.jar 1> /dev/null 2>&1 && ls ionquant/IonQuant*.jar 1> /dev/null 2>&1 && ls diatracer/diaTracer*.jar 1> /dev/null 2>&1; then
        TOOLS_EXIST=true
        echo "\${GREEN}=========================================================================="
        echo "Tools already downloaded and cached!"
        echo "MSFragger: \$(ls msfragger/MSFragger*.jar)"
        echo "IonQuant: \$(ls ionquant/IonQuant*.jar)"
        echo "DiaTracer: \$(ls diatracer/diaTracer*.jar)"
        echo "==========================================================================\${RESET}"
        echo ""
        echo "Skipping download. Using cached tools."
    fi

    if [ "\$TOOLS_EXIST" = "false" ]; then
        echo "Tools not found. Download required."
        echo ""

        #
        # STEP 2: Determine mode (interactive vs token-only vs full non-interactive)
        #
        INTERACTIVE="${interactive_mode}"
        TOKEN_ONLY="${token_only_mode}"
        FULL_MODE="${full_mode}"

        if [ "\$TOKEN_ONLY" = "true" ]; then
            #
            # TOKEN-ONLY MODE
            #
            echo "\${CYAN}=========================================================================="
            echo "                    TOKEN-ONLY MODE"
            echo "==========================================================================\${RESET}"
            echo ""
            echo "Using token provided via --fragpipe_token parameter"
            echo "Skipping registration (assuming you already obtained token via curl)"
            echo ""
            USER_TOKEN="${token}"
        elif [ "\$INTERACTIVE" = "true" ]; then
            #
            # INTERACTIVE MODE
            #
            echo "\${YELLOW}=========================================================================="
            echo "                    INTERACTIVE SETUP MODE"
            echo "==========================================================================\${RESET}"
            echo ""
            echo "This wizard will guide you through FragPipe setup."
            echo "You will need to:"
            echo "  1. Accept the license agreement"
            echo "  2. Provide your contact information"
            echo "  3. Receive and enter a verification code from your email"
            echo ""

            # Display license information
            echo "\${CYAN}=========================================================================="
            echo "                    FRAGPIPE LICENSE INFORMATION"
            echo "==========================================================================\${RESET}"
            cat ${projectDir}/assets/licensing_information.txt
            echo ""

            # Request license acceptance
            while true; do
                read -p "\${YELLOW}Do you accept these license terms? (yes/no): \${RESET}" LICENSE_RESPONSE
                case "\$LICENSE_RESPONSE" in
                    [Yy]|[Yy][Ee][Ss])
                        echo "\${GREEN}✓ License accepted\${RESET}"
                        break
                        ;;
                    [Nn]|[Nn][Oo])
                        echo "\${RED}ERROR: License not accepted. Cannot proceed.\${RESET}"
                        exit 1
                        ;;
                    *)
                        echo "\${RED}Invalid response. Please enter 'yes' or 'no'.\${RESET}"
                        ;;
                esac
            done

            # Display additional license information
            echo ""
            echo "\${CYAN}=========================================================================="
            echo "           ADDITIONAL MSFRAGGER LICENSE INFORMATION"
            echo "==========================================================================\${RESET}"
            cat ${projectDir}/assets/msfragger_licensing_information.txt
            echo ""

            # Request additional license acceptance
            while true; do
                read -p "\${YELLOW}Do you accept these additional license terms? (yes/no): \${RESET}" LICENSE_RESPONSE2
                case "\$LICENSE_RESPONSE2" in
                    [Yy]|[Yy][Ee][Ss])
                        echo "\${GREEN}✓ Additional license accepted\${RESET}"
                        break
                        ;;
                    [Nn]|[Nn][Oo])
                        echo "\${RED}ERROR: Additional license not accepted. Cannot proceed.\${RESET}"
                        exit 1
                        ;;
                    *)
                        echo "\${RED}Invalid response. Please enter 'yes' or 'no'.\${RESET}"
                        ;;
                esac
            done

            echo ""
            echo "\${CYAN}=========================================================================="
            echo "                    USER INFORMATION"
            echo "==========================================================================\${RESET}"
            echo ""
            echo "Please enter your contact information for registration:"
            echo ""
            read -p "First Name: " USER_FIRST_NAME
            read -p "Last Name: " USER_LAST_NAME
            read -p "Email: " USER_EMAIL
            read -p "Institution/Organization: " USER_INSTITUTION
            echo ""
            echo "Registration information:"
            echo "  Name: \$USER_FIRST_NAME \$USER_LAST_NAME"
            echo "  Email: \$USER_EMAIL"
            echo "  Institution: \$USER_INSTITUTION"
            echo ""

        else
            #
            # NON-INTERACTIVE MODE
            #
            echo "\${CYAN}=========================================================================="
            echo "                    NON-INTERACTIVE SETUP MODE"
            echo "==========================================================================\${RESET}"
            echo ""
            echo "Using parameters provided via command line:"
            echo "  Name: ${first_name} ${last_name}"
            echo "  Email: ${email}"
            echo "  Institution: ${institution}"
            echo "  License Accepted: ${license_accept}"
            echo ""

            USER_FIRST_NAME="${first_name}"
            USER_LAST_NAME="${last_name}"
            USER_EMAIL="${email}"
            USER_INSTITUTION="${institution}"
        fi

        #
        # STEP 3: Get latest version and send registration (skip registration for token-only mode)
        #
        # Get latest MSFragger version (needed for all modes for download URLs)
        echo "\${YELLOW}Querying latest MSFragger version...\${RESET}"
        MSFRAGGER_VERSION=\$(curl -s https://msfragger-upgrader.nesvilab.org/upgrader/latest_version.php)
        echo "Latest MSFragger version: \$MSFRAGGER_VERSION"
        echo ""

        # Send registration to Nesvilab (skip for token-only mode since user already did this)
        if [ "\$TOKEN_ONLY" = "false" ]; then
            echo "\${YELLOW}Registering with Nesvilab upgrader server...\${RESET}"

            curl --location --request POST \\
                'https://msfragger-upgrader.nesvilab.org/upgrader/upgrade_download.php' \\
                --form 'transfer="academic"' \\
                --form 'agreement2="true"' \\
                --form 'agreement3="true"' \\
                --form "first_name=\$USER_FIRST_NAME" \\
                --form "last_name=\$USER_LAST_NAME" \\
                --form "email=\$USER_EMAIL" \\
                --form "organization=\$USER_INSTITUTION" \\
                --form "download=\${MSFRAGGER_VERSION}\\\$zip" \\
                --form 'is_fragpipe="true"' \\
                > /dev/null 2>&1

            echo "\${GREEN}✓ Registration sent successfully!\${RESET}"
            echo ""
        fi

        #
        # STEP 4: Get authentication token (skip for token-only mode - already set)
        #
        if [ "\$TOKEN_ONLY" = "false" ]; then
            if [ "\$INTERACTIVE" = "true" ]; then
                echo "\${YELLOW}=========================================================================="
                echo "                    AUTHENTICATION TOKEN REQUIRED"
                echo "==========================================================================\${RESET}"
                echo ""
                echo "A 6-digit verification code has been sent to: \$USER_EMAIL"
                echo "\${CYAN}Please check your email (including spam folder)\${RESET}"
                echo ""

                # Prompt for token with validation
                while true; do
                    read -p "Enter the 6-digit authentication code: " USER_TOKEN
                    if [[ "\$USER_TOKEN" =~ ^[0-9]{6}\$ ]]; then
                        echo "\${GREEN}✓ Token accepted: \$USER_TOKEN\${RESET}"
                        break
                    else
                        echo "\${RED}✗ Invalid token. Please enter exactly 6 digits.\${RESET}"
                    fi
                done
            else
                echo "Using token from parameters: ${token}"
                USER_TOKEN="${token}"
            fi
        fi

        echo ""

        #
        # STEP 5: Download tools using token
        #
        echo "\${YELLOW}Downloading MSFragger and IonQuant...\${RESET}"
        echo ""
        echo "Downloading MSFragger..."
        MSFRAGGER_VERSION_ENCODED=\${MSFRAGGER_VERSION// /%20}

        curl --fail --silent --show-error --location --output msfragger.zip \\
            "https://msfragger-upgrader.nesvilab.org/upgrader/download.php?token=\$USER_TOKEN&download=\${MSFRAGGER_VERSION_ENCODED}%24zip"

        if [ ! -f msfragger.zip ] || [ ! -s msfragger.zip ]; then
            echo "\${RED}ERROR: Failed to download MSFragger\${RESET}"
            echo "Please check your token is correct: \$USER_TOKEN"
            echo "Token should be sent to: \$USER_EMAIL"
            exit 1
        fi

        # Check if download worked (downloaded software is in ZIP format); if not prompt for new token since expired
        if ! file msfragger.zip | grep -q "Zip archive"; then
            echo ""
            echo -e "\${RED}════════════════════════════════════════════════════════════════════════\${RESET}"
            echo -e "\${RED}                    AUTHENTICATION TOKEN EXPIRED\${RESET}"
            echo -e "\${RED}════════════════════════════════════════════════════════════════════════\${RESET}"
            echo ""
            echo -e "\${YELLOW}Your 6-digit verification code (\$USER_TOKEN) has expired or is invalid.\${RESET}"
            echo ""
            echo "To fix this issue:"
            echo ""

            if [ "\$TOKEN_ONLY" = "true" ]; then
                echo -e "  \${CYAN}1.\${RESET} Request a new token by running this curl command:"
                echo ""
                echo -e "\${GREEN}"
                echo "     curl --location --request POST \\\\"
                echo "       'https://msfragger-upgrader.nesvilab.org/upgrader/upgrade_download.php' \\\\"
                echo "       --form 'transfer=\"academic\"' \\\\"
                echo "       --form 'agreement2=\"true\"' \\\\"
                echo "       --form 'agreement3=\"true\"' \\\\"
                echo "       --form 'first_name=YOUR_FIRST_NAME' \\\\"
                echo "       --form 'last_name=YOUR_LAST_NAME' \\\\"
                echo "       --form 'email=YOUR_EMAIL' \\\\"
                echo "       --form 'organization=YOUR_INSTITUTION' \\\\"
                echo "       --form 'download=4.4.1\\\\\\\$zip' \\\\"
                echo "       --form 'is_fragpipe=\"true\"'"
                echo -e "\${RESET}"
                echo ""
                echo -e "  \${CYAN}2.\${RESET} Check your email for the new 6-digit token"
                echo ""
                echo -e "  \${CYAN}3.\${RESET} Update your pipeline command with the new token:"
                echo ""
                echo -e "     \${GREEN}--fragpipe_token XXXXXX\${RESET}"
            else
                echo -e "  \${CYAN}1.\${RESET} Check your email (\$USER_EMAIL) for a \${CYAN}NEW\${RESET} verification code"
                echo -e "     \${YELLOW}(The code expires quickly - usually within minutes to hours)\${RESET}"
                echo ""
                echo -e "  \${CYAN}2.\${RESET} If you don't have a recent email, request a new token by running:"
                echo ""
                echo -e "\${GREEN}"
                echo "     curl --location --request POST \\\\"
                echo "       'https://msfragger-upgrader.nesvilab.org/upgrader/upgrade_download.php' \\\\"
                echo "       --form 'transfer=\"academic\"' \\\\"
                echo "       --form 'agreement2=\"true\"' \\\\"
                echo "       --form 'agreement3=\"true\"' \\\\"
                echo "       --form 'first_name=\$USER_FIRST_NAME' \\\\"
                echo "       --form 'last_name=\$USER_LAST_NAME' \\\\"
                echo "       --form 'email=\$USER_EMAIL' \\\\"
                echo "       --form 'organization=\$USER_INSTITUTION' \\\\"
                echo "       --form 'download=\${MSFRAGGER_VERSION}\\\$zip' \\\\"
                echo "       --form 'is_fragpipe=\"true\"'"
                echo -e "\${RESET}"
                echo ""
                echo -e "  \${CYAN}3.\${RESET} Update your pipeline command with the new token:"
                echo ""
                echo -e "     \${GREEN}--fragpipe_token XXXXXX\${RESET}"
            fi
            echo ""
            echo -e "  \${CYAN}4.\${RESET} Resume the pipeline run:"
            echo ""
            echo -e "     \${GREEN}nextflow run <your_pipeline> -resume [other parameters]\${RESET}"
            echo ""
            echo -e "\${RED}════════════════════════════════════════════════════════════════════════\${RESET}"
            echo ""
            exit 1
        fi

        unzip -q msfragger.zip
        find . -maxdepth 3 -name "MSFragger*.jar" -exec mv {} msfragger/ \\;
        rm -f msfragger.zip

        if [ ! -f msfragger/MSFragger*.jar ]; then
            echo "\${RED}ERROR: MSFragger JAR not found after extraction\${RESET}"
            ls -R
            exit 1
        fi

        echo "\${GREEN}✓ MSFragger downloaded: \$(ls msfragger/MSFragger*.jar)\${RESET}"

        echo "Downloading IonQuant..."
        IONQUANT_VERSION=\$(curl -s https://msfragger-upgrader.nesvilab.org/ionquant/latest_version.php)
        curl --fail --silent --show-error --location --output ionquant.zip \\
            "https://msfragger-upgrader.nesvilab.org/ionquant/download.php?token=\$USER_TOKEN&download=\${IONQUANT_VERSION}%24zip"

        if [ ! -f ionquant.zip ] || [ ! -s ionquant.zip ]; then
            echo "\${RED}ERROR: Failed to download IonQuant\${RESET}"
            echo "Please check your token is correct: \$USER_TOKEN"
            exit 1
        fi

        # Check if token worked
        if ! file ionquant.zip | grep -q "Zip archive"; then
            echo ""
            echo -e "\${RED}════════════════════════════════════════════════════════════════════════\${RESET}"
            echo -e "\${RED}                    AUTHENTICATION TOKEN EXPIRED\${RESET}"
            echo -e "\${RED}════════════════════════════════════════════════════════════════════════\${RESET}"
            echo ""
            echo -e "\${YELLOW}Your 6-digit verification code (\$USER_TOKEN) has expired or is invalid.\${RESET}"
            echo ""
            echo "The token expired while downloading IonQuant."
            echo "Please obtain a fresh verification code and try again."
            echo ""
            echo "See instructions above for how to get a new token."
            echo ""
            echo -e "\${RED}════════════════════════════════════════════════════════════════════════\${RESET}"
            echo ""
            exit 1
        fi

        unzip -q ionquant.zip
        find . -maxdepth 3 -name "IonQuant*.jar" -exec mv {} ionquant/ \\;
        rm -f ionquant.zip

        if [ ! -f ionquant/IonQuant*.jar ]; then
            echo "\${RED}ERROR: IonQuant JAR not found after extraction\${RESET}"
            ls -R
            exit 1
        fi

        echo "\${GREEN}✓ IonQuant downloaded: \$(ls ionquant/IonQuant*.jar)\${RESET}"

        echo "Downloading DiaTracer..."
        DIATRACER_VERSION=\$(curl -s https://msfragger-upgrader.nesvilab.org/diatracer/latest_version.php)
        curl --fail --silent --show-error --location --output diatracer.zip \
            "https://msfragger-upgrader.nesvilab.org/diatracer/download.php?token=\$USER_TOKEN&download=\${DIATRACER_VERSION}%24zip"

        if [ ! -f diatracer.zip ] || [ ! -s diatracer.zip ]; then
            echo "\${RED}ERROR: Failed to download DiaTracer\${RESET}"
            echo "Please check your token is correct: \$USER_TOKEN"
            exit 1
        fi

        # Check if token worked
        if ! file diatracer.zip | grep -q "Zip archive"; then
            echo ""
            echo -e "\${RED}════════════════════════════════════════════════════════════════════════\${RESET}"
            echo -e "\${RED}                    AUTHENTICATION TOKEN EXPIRED\${RESET}"
            echo -e "\${RED}════════════════════════════════════════════════════════════════════════\${RESET}"
            echo ""
            echo -e "\${YELLOW}Your 6-digit verification code (\$USER_TOKEN) has expired or is invalid.\${RESET}"
            echo ""
            echo "The token expired while downloading DiaTracer."
            echo "Please obtain a fresh verification code and try again."
            echo ""
            echo "See instructions above for how to get a new token."
            echo ""
            echo -e "\${RED}════════════════════════════════════════════════════════════════════════\${RESET}"
            echo ""
            exit 1
        fi

        unzip -q diatracer.zip
        find . -maxdepth 3 -name "diaTracer*.jar" -exec mv {} diatracer/ \\;
        rm -f diatracer.zip

        if [ ! -f diatracer/diaTracer*.jar ]; then
            echo "\${RED}ERROR: DiaTracer JAR not found after extraction\${RESET}"
            ls -R
            exit 1
        fi

        echo "\${GREEN}✓ DiaTracer downloaded: \$(ls diatracer/diaTracer*.jar)\${RESET}"
        echo ""
        echo "\${GREEN}=========================================================================="
        echo "FragPipe tools downloaded successfully!"
        echo "==========================================================================\${RESET}"
    fi

    #
    # FINAL VERIFICATION
    #
    echo ""
    echo "Verifying tools..."
    MSFRAGGER_JAR=\$(ls msfragger/MSFragger*.jar)
    IONQUANT_JAR=\$(ls ionquant/IonQuant*.jar)
    DIATRACER_JAR=\$(ls diatracer/diaTracer*.jar)

    if [ ! -f "\$MSFRAGGER_JAR" ] || [ ! -s "\$MSFRAGGER_JAR" ]; then
        echo "\${RED}ERROR: MSFragger JAR is missing or empty\${RESET}"
        exit 1
    fi

    if [ ! -f "\$IONQUANT_JAR" ] || [ ! -s "\$IONQUANT_JAR" ]; then
        echo "\${RED}ERROR: IonQuant JAR is missing or empty\${RESET}"
        exit 1
    fi

    if [ ! -f "\$DIATRACER_JAR" ] || [ ! -s "\$DIATRACER_JAR" ]; then
        echo "\${RED}ERROR: DiaTracer JAR is missing or empty\${RESET}"
        exit 1
    fi

    echo "\${GREEN}✓ All tools verified and ready!\${RESET}"
    echo ""
    echo "Tool locations (cached for future runs):"
    echo "  MSFragger: \$MSFRAGGER_JAR"
    echo "  IonQuant: \$IONQUANT_JAR"
    echo "  DiaTracer: \$DIATRACER_JAR"
    echo ""
    echo "These tools will be reused in future pipeline runs."
    echo "\${GREEN}=========================================================================="
    echo "Setup complete!"
    echo "==========================================================================\${RESET}"
    echo ""

    # Generate versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curl: \$(curl --version 2>&1 | head -n 1 | sed 's/^curl //; s/ .*//')
        msfragger: \$(basename \$MSFRAGGER_JAR | grep -oP 'MSFragger-\\K[0-9.]+' || echo "unknown")
        ionquant: \$(basename \$IONQUANT_JAR | grep -oP 'IonQuant-\\K[0-9.]+' || echo "unknown")
        diatracer: \$(basename \$DIATRACER_JAR | grep -oP 'diaTracer-\\K[0-9.]+' || echo "unknown")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p msfragger ionquant diatracer
    touch msfragger/MSFragger-4.1.jar
    touch ionquant/IonQuant-1.10.12.jar
    touch diatracer/diaTracer-1.1.3.jar

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        curl: 8.1.2
        msfragger: 4.1
        ionquant: 1.10.12
        diatracer: 1.1.3
    END_VERSIONS
    """
}
