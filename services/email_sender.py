"""
Send email via Microsoft Graph API using OAuth2 client credentials.

Required env vars:
    AZURE_TENANT_ID
    AZURE_CLIENT_ID
    AZURE_CLIENT_SECRET
    EMAIL_FROM          (the mailbox to send from, e.g. DONOTREPLY@genodynlabs.com)
"""
import os
import base64
import msal
import requests


GRAPH_BASE = "https://graph.microsoft.com/v1.0"
SCOPE = ["https://graph.microsoft.com/.default"]


def _get_access_token() -> str:
    """Obtain an app-only access token from Microsoft Entra (Azure AD)."""
    tenant_id     = os.environ["AZURE_TENANT_ID"]
    client_id     = os.environ["AZURE_CLIENT_ID"]
    client_secret = os.environ["AZURE_CLIENT_SECRET"]

    authority = f"https://login.microsoftonline.com/{tenant_id}"
    app = msal.ConfidentialClientApplication(
        client_id=client_id,
        client_credential=client_secret,
        authority=authority,
    )

    result = app.acquire_token_for_client(scopes=SCOPE)
    if "access_token" not in result:
        raise RuntimeError(
            f"Failed to get Azure access token: "
            f"{result.get('error')} - {result.get('error_description')}"
        )
    return result["access_token"]


def send_report_email(
    to_address: str,
    user_display_name: str,
    pdf_bytes: bytes,
    pdf_filename: str = "genodynlabs_report.pdf",
) -> None:
    """
    Send the PDF report to the user via Microsoft Graph (sendMail endpoint).
    Raises on failure.
    """
    from_address = os.environ["EMAIL_FROM"]
    token = _get_access_token()

    greeting_name = user_display_name or "there"

    body_html = f"""
    <p>Hi {greeting_name},</p>
    <p>Your GenoDynLabs genomic analysis report is attached as a PDF.</p>
    <p>This report contains your annotated SNP variants identified from your
       uploaded raw DNA file, cross-referenced with the SNPedia knowledge base.</p>
    <p>If you have any questions about your results, just reply to this email.</p>
    <p>— The GenoDynLabs team</p>
    """

    payload = {
        "message": {
            "subject": "Your GenoDynLabs Genomic Report",
            "body": {
                "contentType": "HTML",
                "content": body_html,
            },
            "toRecipients": [
                {"emailAddress": {"address": to_address}}
            ],
            "attachments": [
                {
                    "@odata.type": "#microsoft.graph.fileAttachment",
                    "name": pdf_filename,
                    "contentType": "application/pdf",
                    "contentBytes": base64.b64encode(pdf_bytes).decode("ascii"),
                }
            ],
        },
        "saveToSentItems": "false",
    }

    url = f"{GRAPH_BASE}/users/{from_address}/sendMail"
    headers = {
        "Authorization": f"Bearer {token}",
        "Content-Type": "application/json",
    }

    resp = requests.post(url, headers=headers, json=payload, timeout=30)
    if resp.status_code not in (200, 202):
        raise RuntimeError(
            f"Graph sendMail failed ({resp.status_code}): {resp.text[:500]}"
        )