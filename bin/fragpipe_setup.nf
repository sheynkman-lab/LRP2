#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * FragPipe Interactive Setup Functions
 * Based on FragNFlow's config_tools_init.nf
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

def displayLicensingInformation() {
    println ""
    println "=========================================================================="
    println "                    FRAGPIPE LICENSE INFORMATION"
    println "=========================================================================="
    println file(projectDir + '/assets/licensing_information.txt').text
}

def displayAdditionalLicensingInformation() {
    println ""
    println "=========================================================================="
    println "           ADDITIONAL MSFRAGGER LICENSE INFORMATION"
    println "=========================================================================="
    println file(projectDir + '/assets/msfragger_licensing_information.txt').text
}

def checkToolsAvailability() {
    def tools_dir = params.fragpipe_tools_dir
    def msfragger_exists = file("${tools_dir}/msfragger").exists() &&
                          file("${tools_dir}/msfragger").isDirectory() &&
                          file("${tools_dir}/msfragger").list().findAll { it.endsWith('.jar') }.size() > 0
    def ionquant_exists = file("${tools_dir}/ionquant").exists() &&
                         file("${tools_dir}/ionquant").isDirectory() &&
                         file("${tools_dir}/ionquant").list().findAll { it.endsWith('.jar') }.size() > 0

    println ""
    println "=========================================================================="
    println "                       TOOLS AVAILABILITY CHECK"
    println "=========================================================================="
    println "Tools directory: ${tools_dir}"
    println ""
    println msfragger_exists ? "${GREEN}✅ MSFragger : available${RESET}" : "${RED}❌ MSFragger : not available${RESET}"
    println ionquant_exists ? "${GREEN}✅ IonQuant : available${RESET}" : "${RED}❌ IonQuant : not available${RESET}"
    println ""

    return [
        msfragger: msfragger_exists,
        ionquant: ionquant_exists,
        all_available: msfragger_exists && ionquant_exists
    ]
}

def collectUserInformation() {
    println ""
    println "${YELLOW}=========================================================================="
    println "                    FRAGPIPE SETUP - USER INFORMATION"
    println "==========================================================================${RESET}"
    println ""
    println "Tools are not available. You need to download them from Nesvilab."
    println "This requires registration and accepting the license agreement."
    println ""
    println "${CYAN}Please enter your contact information:${RESET}"
    println ""

    print "First Name: "
    def first_name = System.in.newReader().readLine()

    print "Last Name: "
    def last_name = System.in.newReader().readLine()

    print "Email: "
    def email = System.in.newReader().readLine()

    print "Institution/Organization: "
    def institution = System.in.newReader().readLine()

    println ""

    return [
        first_name: first_name,
        last_name: last_name,
        email: email,
        institution: institution
    ]
}

def requestLicenseAcceptance() {
    displayLicensingInformation()

    println ""
    print "${YELLOW}Do you accept these license terms? (yes/no): ${RESET}"
    def response = System.in.newReader().readLine()

    if (response.toLowerCase() ==~ /yes|y/) {
        println "${GREEN}✓ License accepted${RESET}"
        return true
    } else {
        error "${RED}License not accepted. Cannot proceed with FragPipe setup.${RESET}"
    }
}

def requestAdditionalLicenseAcceptance() {
    displayAdditionalLicensingInformation()

    println ""
    print "${YELLOW}Do you accept these additional license terms? (yes/no): ${RESET}"
    def response = System.in.newReader().readLine()

    if (response.toLowerCase() ==~ /yes|y/) {
        println "${GREEN}✓ Additional license accepted${RESET}"
        return true
    } else {
        error "${RED}Additional license not accepted. Cannot proceed with FragPipe setup.${RESET}"
    }
}

def promptForToken() {
    println ""
    println "${YELLOW}=========================================================================="
    println "                    AUTHENTICATION TOKEN REQUIRED"
    println "==========================================================================${RESET}"
    println ""
    println "A 6-digit verification code has been sent to your email."
    println "${CYAN}Please check your email (including spam folder)${RESET}"
    println ""

    def token
    while (true) {
        print "Enter the 6-digit authentication code: "
        token = System.in.newReader().readLine()

        if (token ==~ /\d{6}/) {
            println "${GREEN}✓ Token accepted: ${token}${RESET}"
            break
        } else {
            println "${RED}✗ Invalid input.${RESET} Please enter a 6-digit code."
        }
    }

    return token
}

def runInteractiveSetup() {
    println ""
    println "${CYAN}=========================================================================="
    println "                    FRAGPIPE INTERACTIVE SETUP"
    println "==========================================================================${RESET}"
    println ""
    println "This wizard will guide you through setting up FragPipe."
    println "You will need to:"
    println "  1. Accept the license agreement"
    println "  2. Provide your contact information"
    println "  3. Receive and enter a verification code from your email"
    println ""

    // Check if tools already exist
    def tools_status = checkToolsAvailability()

    if (tools_status.all_available) {
        println "${GREEN}=========================================================================="
        println "All tools are already downloaded and available!"
        println "==========================================================================${RESET}"
        println ""
        println "You can skip the setup and run the pipeline directly."
        println ""
        return [
            skip_download: true,
            tools_ready: true
        ]
    }

    // Collect user information
    def user_info = collectUserInformation()

    // Request license acceptance
    def license_accepted = requestLicenseAcceptance()
    def additional_license_accepted = requestAdditionalLicenseAcceptance()

    println ""
    println "${GREEN}=========================================================================="
    println "Registration information will be sent to Nesvilab..."
    println "==========================================================================${RESET}"
    println ""
    println "Name: ${user_info.first_name} ${user_info.last_name}"
    println "Email: ${user_info.email}"
    println "Institution: ${user_info.institution}"
    println ""

    return [
        skip_download: false,
        first_name: user_info.first_name,
        last_name: user_info.last_name,
        email: user_info.email,
        institution: user_info.institution,
        license_accept: license_accepted
    ]
}
