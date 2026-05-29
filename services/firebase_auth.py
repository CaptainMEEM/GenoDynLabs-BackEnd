"""
Firebase Admin SDK initialization and token verification.

Reads service account JSON from FIREBASE_CREDENTIALS_JSON env var
(either as raw JSON string or as a path to a file).
"""
import os
import json
import firebase_admin
from firebase_admin import credentials, auth as fb_auth
from functools import wraps
from flask import request, jsonify


_initialized = False


def init_firebase():
    """Initialize the Firebase Admin SDK once."""
    global _initialized
    if _initialized:
        return

    raw = os.environ.get("FIREBASE_CREDENTIALS_JSON")
    if not raw:
        raise RuntimeError(
            "FIREBASE_CREDENTIALS_JSON env var is not set. "
            "Paste the contents of your serviceAccountKey.json into Railway's env vars."
        )

    # Accept either raw JSON or a path to a file
    if raw.strip().startswith("{"):
        cred_dict = json.loads(raw)
        cred = credentials.Certificate(cred_dict)
    else:
        cred = credentials.Certificate(raw)

    firebase_admin.initialize_app(cred)
    _initialized = True


def verify_token(id_token: str) -> dict:
    """
    Verify a Firebase ID token. Returns the decoded token (with uid, email, etc.)
    or raises an exception if invalid.
    """
    if not _initialized:
        init_firebase()
    return fb_auth.verify_id_token(id_token)


def require_auth(f):
    """
    Flask decorator. Rejects the request with 401 if there's no valid
    Firebase ID token in the Authorization header. On success, injects
    the decoded token as the first argument.

    Usage:
        @app.route("/secret")
        @require_auth
        def secret(user):
            return f"hello {user['email']}"
    """
    @wraps(f)
    def wrapper(*args, **kwargs):
        header = request.headers.get("Authorization", "")
        if not header.startswith("Bearer "):
            return jsonify({"error": "Missing or invalid Authorization header"}), 401

        token = header.split(" ", 1)[1].strip()
        try:
            user = verify_token(token)
        except Exception as e:
            return jsonify({"error": f"Invalid token: {e}"}), 401

        return f(user, *args, **kwargs)
    return wrapper
