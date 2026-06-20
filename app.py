"""
PubMed 논문 검색 웹 앱 (Streamlit)
실행: streamlit run app.py
"""

import time
import io
import re
import xml.etree.ElementTree as ET
import streamlit as st
import pandas as pd
from Bio import Entrez, Medline

# ─────────────────────────────────────────────
# NCBI 접속 정보
# ─────────────────────────────────────────────
Entrez.email   = "sbhong2008@gmail.com"
Entrez.api_key = "09766502d4a2b055f09426a9658e5bf33808"

EMAIL_RE = re.compile(r"[\w.\-+]+@[\w.\-]+\.[a-zA-Z]{2,}")


# ─────────────────────────────────────────────
# 검색 관련 함수
# ─────────────────────────────────────────────
def _wildcard_terms(text: str, field: str) -> str:
    words = text.strip().split()
    return " AND ".join(f"{w}*[{field}]" for w in words)


def build_query(keyword: str = "", author: str = "", title: str = "") -> str:
    parts = []
    if keyword.strip():
        words = keyword.strip().split()
        parts.append(" AND ".join(f"{w}*" for w in words))
    if author.strip():
        words = author.strip().split()
        parts.append(" AND ".join(f"{w}*[Author]" for w in words))
    if title.strip():
        parts.append(_wildcard_terms(title, "Title"))
    if not parts:
        raise ValueError("키워드, 저자, 제목 중 하나 이상 입력하세요.")
    return " AND ".join(parts)


def search_pubmed(query: str, retmax: int):
    handle = Entrez.esearch(db="pubmed", term=query, retmax=retmax)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"], int(record["Count"])


def fetch_details(id_list: list[str], progress_bar, batch_size: int = 10) -> list[dict]:
    records = []
    total = len(id_list)
    for start in range(0, total, batch_size):
        batch_ids = id_list[start:start + batch_size]
        handle = Entrez.efetch(
            db="pubmed", id=",".join(batch_ids),
            rettype="medline", retmode="text"
        )
        records.extend(list(Medline.parse(handle)))
        handle.close()
        progress_bar.progress(min((start + batch_size) / total, 1.0))
        time.sleep(0.4 if Entrez.api_key else 1.0)
    return records


def _all_words_match(text: str, words: list[str]) -> bool:
    t = text.lower()
    return all(any(tok.startswith(w) for tok in t.split()) for w in words)


def filter_records(records, author_filter: str, title_filter: str) -> list[dict]:
    af_words = [w.lower() for w in author_filter.split() if w]
    tf_words = [w.lower() for w in title_filter.split() if w]
    if not af_words and not tf_words:
        return records
    filtered = []
    for rec in records:
        if af_words and not _all_words_match(" ".join(rec.get("AU", [])), af_words):
            continue
        if tf_words and not _all_words_match(rec.get("TI", ""), tf_words):
            continue
        filtered.append(rec)
    return filtered


def parse_records(records: list[dict]) -> pd.DataFrame:
    rows = []
    for rec in records:
        authors = "; ".join(rec.get("AU", []))
        doi = ""
        for aid in rec.get("AID", []):
            if aid.endswith("[doi]"):
                doi = f"https://doi.org/{aid.replace('[doi]', '').strip()}"
                break
        rows.append({
            "PMID":   rec.get("PMID", ""),
            "제목":   rec.get("TI",   ""),
            "저자":   authors,
            "출판일": rec.get("DP",   ""),
            "초록":   rec.get("AB",   ""),
            "DOI":    doi,
        })
    return pd.DataFrame(rows, columns=["PMID", "제목", "저자", "출판일", "초록", "DOI"])


def to_csv_bytes(df: pd.DataFrame) -> bytes:
    buf = io.BytesIO()
    df.to_csv(buf, index=False, encoding="utf-8-sig")
    return buf.getvalue()


# ─────────────────────────────────────────────
# 교신저자 조회 함수
# ─────────────────────────────────────────────
def doi_to_pmid(doi: str) -> str | None:
    """DOI → PMID 변환. URL 접두어는 자동 제거합니다."""
    doi = re.sub(r"https?://doi\.org/", "", doi).strip()
    handle = Entrez.esearch(db="pubmed", term=f"{doi}[DOI]")
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"][0] if record["IdList"] else None


def fetch_corresponding_author(pmid: str) -> dict:
    """
    PubMed XML에서 교신저자 정보를 추출합니다.

    추출 전략 (우선순위 순):
    1. 소속 문자열에 'correspond' / 'electronic address' 키워드 포함 저자
    2. 소속 문자열에 이메일이 포함된 저자
    3. 마지막 저자 (관례상 교신저자인 경우가 많음)
    """
    handle = Entrez.efetch(db="pubmed", id=pmid, rettype="xml", retmode="xml")
    xml_bytes = handle.read()
    handle.close()
    time.sleep(0.4 if Entrez.api_key else 1.0)

    root = ET.fromstring(xml_bytes)

    # 논문 기본 정보
    title   = root.findtext(".//ArticleTitle") or ""
    pub_date_node = root.find(".//PubDate")
    pub_date = ""
    if pub_date_node is not None:
        parts = [pub_date_node.findtext(t) or "" for t in ("Year", "Month", "Day")]
        pub_date = " ".join(p for p in parts if p)

    # 저자 목록 파싱
    authors_data = []
    for au in root.findall(".//AuthorList/Author[@ValidYN='Y']"):
        last  = au.findtext("LastName")  or ""
        fore  = au.findtext("ForeName")  or ""
        init  = au.findtext("Initials")  or ""
        name  = f"{last} {fore}".strip() or f"{last} {init}".strip()

        affils = [
            aff.text.strip()
            for aff in au.findall(".//AffiliationInfo/Affiliation")
            if aff.text
        ]
        emails = []
        for aff in affils:
            emails.extend(EMAIL_RE.findall(aff))
        # <Identifier Source="ORCID"> 등 추가 식별자
        orcid = ""
        for ident in au.findall(".//Identifier"):
            if ident.get("Source", "").upper() == "ORCID":
                orcid = (ident.text or "").replace("http://orcid.org/", "").strip()

        authors_data.append({
            "name":    name,
            "affils":  affils,
            "emails":  emails,
            "orcid":   orcid,
        })

    # ── 교신저자 판별 ──────────────────────────
    corresponding = None

    # 전략 1: 소속에 'correspond' 또는 'electronic address' 포함
    for au in authors_data:
        for aff in au["affils"]:
            if re.search(r"correspond|electronic\s+address", aff, re.I):
                corresponding = au
                break
        if corresponding:
            break

    # 전략 2: 소속에 이메일이 포함된 첫 번째 저자
    if not corresponding:
        for au in authors_data:
            if au["emails"]:
                corresponding = au
                break

    # 전략 3: 마지막 저자 (관례)
    if not corresponding and authors_data:
        corresponding = authors_data[-1]
        corresponding["_fallback"] = True

    return {
        "pmid":          pmid,
        "title":         title,
        "pub_date":      pub_date,
        "corresponding": corresponding,
        "all_authors":   authors_data,
    }


# ─────────────────────────────────────────────
# Streamlit UI
# ─────────────────────────────────────────────
st.set_page_config(page_title="PubMed 논문 검색기", page_icon="🔬", layout="wide")
st.title("🔬 PubMed 논문 검색기")
st.caption("키워드 · 저자 · 제목 검색 및 교신저자 조회")

# ══════════════════════════════════════════════
# 사이드바 — 검색 조건
# ══════════════════════════════════════════════
with st.sidebar:
    st.header("논문 검색 조건")
    keyword = st.text_input("키워드 (전체 필드)", placeholder="예: CRISPR cancer therapy")
    author  = st.text_input("저자명", placeholder="예: Kim J  /  Zhang Yong")
    title   = st.text_input("논문 제목 키워드", placeholder="예: gene editing")
    retmax  = st.slider("최대 수집 건수", min_value=5, max_value=200, value=50, step=5)
    search_btn = st.button("🔍 검색", use_container_width=True, type="primary")

    st.divider()
    st.markdown(
        "**검색 팁**\n"
        "- 조건을 여러 개 입력하면 **AND** 로 결합됩니다.\n"
        "- 단어 **일부만** 입력해도 매칭됩니다.\n"
        "  예) `CRISP` → CRISPR 매칭\n"
        "  예) `Kim` → Kim J / Kim JS 매칭\n"
        "- 여러 단어를 입력하면 **모두 포함**된 결과만 출력합니다."
    )

# ══════════════════════════════════════════════
# 메인 탭
# ══════════════════════════════════════════════
main_tab1, main_tab2 = st.tabs(["📚 논문 검색", "👤 교신저자 조회"])


# ── 탭 1: 논문 검색 ───────────────────────────
with main_tab1:
    if search_btn:
        if not any([keyword.strip(), author.strip(), title.strip()]):
            st.warning("키워드, 저자, 제목 중 하나 이상 입력하세요.")
            st.stop()

        try:
            query = build_query(keyword=keyword, author=author, title=title)
        except ValueError as e:
            st.error(str(e))
            st.stop()

        st.info(f"**검색 쿼리:** `{query}`")

        with st.spinner("PubMed에서 검색 중..."):
            try:
                id_list, total_count = search_pubmed(query, retmax)
            except Exception as e:
                st.error(f"검색 오류: {e}")
                st.stop()

        if not id_list:
            st.warning("검색 결과가 없습니다.")
            st.stop()

        st.success(f"총 **{total_count:,}건** 중 **{len(id_list)}건** 수집")

        st.write("논문 상세 정보 수집 중...")
        progress = st.progress(0)
        try:
            raw_records = fetch_details(id_list, progress)
        except Exception as e:
            st.error(f"데이터 수집 오류: {e}")
            st.stop()
        progress.empty()

        filtered = filter_records(raw_records, author_filter=author, title_filter=title)
        if not filtered:
            st.warning("필터 조건에 맞는 논문이 없습니다.")
            st.stop()

        df = parse_records(filtered)

        col1, col2, col3 = st.columns(3)
        col1.metric("수집 논문", f"{len(df)}건")
        col2.metric("저자 수 (중복 포함)", f"{df['저자'].str.split(';').explode().nunique()}명")
        col3.metric("DOI 보유", f"{df['DOI'].astype(bool).sum()}건")

        st.divider()

        sub1, sub2 = st.tabs(["📋 논문 목록", "📄 초록 보기"])

        with sub1:
            st.dataframe(
                df[["PMID", "제목", "저자", "출판일", "DOI"]],
                use_container_width=True,
                column_config={"DOI": st.column_config.LinkColumn("DOI 링크")},
                hide_index=True,
            )

        with sub2:
            for _, row in df.iterrows():
                label = f"[{row['출판일']}]  {row['제목'][:90]}{'...' if len(row['제목']) > 90 else ''}"
                with st.expander(label):
                    st.markdown(f"**저자:** {row['저자']}")
                    st.markdown(f"**DOI:** {row['DOI'] or '없음'}")
                    st.markdown("**초록:**")
                    st.write(row["초록"] or "초록 없음")

        st.divider()
        st.download_button(
            label="⬇️ CSV 다운로드",
            data=to_csv_bytes(df),
            file_name="pubmed_results.csv",
            mime="text/csv",
            use_container_width=True,
        )


# ── 탭 2: 교신저자 조회 ───────────────────────
with main_tab2:
    st.subheader("교신저자 조회")
    st.caption("PMID 또는 DOI를 입력하면 교신저자 이름 · 소속 · 이메일을 추출합니다.")

    lookup_input = st.text_input(
        "PMID 또는 DOI 입력",
        placeholder="예: 38476528  또는  https://doi.org/10.1038/s41586-024-07618-3"
    )
    lookup_btn = st.button("🔎 교신저자 조회", type="primary")

    if lookup_btn:
        raw = lookup_input.strip()
        if not raw:
            st.warning("PMID 또는 DOI를 입력하세요.")
            st.stop()

        # DOI vs PMID 자동 판별
        if re.search(r"10\.\d{4,}", raw):          # DOI 패턴
            with st.spinner("DOI → PMID 변환 중..."):
                pmid = doi_to_pmid(raw)
            if not pmid:
                st.error("해당 DOI로 PubMed에서 논문을 찾을 수 없습니다.")
                st.stop()
            st.info(f"DOI → PMID 변환 완료: **{pmid}**")
        else:
            pmid = re.sub(r"\D", "", raw)          # 숫자만 추출
            if not pmid:
                st.error("올바른 PMID 또는 DOI를 입력하세요.")
                st.stop()

        with st.spinner(f"PMID {pmid} 교신저자 정보 조회 중..."):
            try:
                result = fetch_corresponding_author(pmid)
            except Exception as e:
                st.error(f"조회 오류: {e}")
                st.stop()

        corr = result["corresponding"]

        # ── 논문 기본 정보 ──
        st.markdown(f"### {result['title']}")
        st.caption(f"PMID: {result['pmid']}  |  출판일: {result['pub_date']}")
        st.divider()

        # ── 교신저자 카드 ──
        if corr:
            is_fallback = corr.get("_fallback", False)
            label = "마지막 저자 (교신저자 명시 없음)" if is_fallback else "교신저자"
            st.markdown(f"#### {label}")

            info_col, detail_col = st.columns([1, 2])
            with info_col:
                st.markdown(f"**이름**")
                st.write(corr["name"] or "—")
                st.markdown(f"**이메일**")
                if corr["emails"]:
                    for em in corr["emails"]:
                        st.write(em)
                else:
                    st.write("—")
                if corr["orcid"]:
                    st.markdown("**ORCID**")
                    st.write(corr["orcid"])

            with detail_col:
                st.markdown("**소속**")
                if corr["affils"]:
                    for aff in corr["affils"]:
                        st.write(aff)
                else:
                    st.write("—")
        else:
            st.warning("교신저자 정보를 찾을 수 없습니다.")

        # ── 전체 저자 목록 ──
        with st.expander("전체 저자 목록 보기"):
            for i, au in enumerate(result["all_authors"], 1):
                emails_str = ", ".join(au["emails"]) if au["emails"] else "—"
                orcid_str  = au["orcid"] or "—"
                aff_str    = "\n".join(au["affils"]) if au["affils"] else "—"
                st.markdown(
                    f"**{i}. {au['name']}**  \n"
                    f"소속: {aff_str}  \n"
                    f"이메일: {emails_str}  |  ORCID: {orcid_str}"
                )
                st.divider()
